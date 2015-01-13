#ifndef _SNARKLIB_PAIRING_HPP_
#define _SNARKLIB_PAIRING_HPP_

#include <cstdint>
#include <gmp.h>
#include <istream>
#include <ostream>
#include <vector>
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "Group.hpp"
#include "ProgressCallback.hpp"
#include "WindowExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Paired group knowledge commitment
//

template <typename GA, typename GB>
class Pairing
{
public:
    Pairing() = default;

    Pairing(const GA& a, const GB& b)
        : m_G(a), m_H(b)
    {}

    // breaks encapsulation
    const GA& G() const { return m_G; }
    const GB& H() const { return m_H; }

    bool operator== (const Pairing& other) const {
        return m_G == other.m_G && m_H == other.m_H;
    }

    bool operator!= (const Pairing& other) const {
        return ! (*this == other);
    }

    bool isZero() const {
        return *this == zero();
    }

    static Pairing zero() {
        return Pairing<GA, GB>(GA::zero(), GB::zero());
    }

    void marshal_out(std::ostream& os) const {
        G().marshal_out(os);
        H().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_G.marshal_in(is) &&
            m_H.marshal_in(is);
    }

private:
    GA m_G;
    GB m_H;
};

////////////////////////////////////////////////////////////////////////////////
// Operator functions
//

template <typename GA, typename GB>
Pairing<GA, GB> operator+ (const Pairing<GA, GB>& a, const Pairing<GA, GB>& b) {
    return Pairing<GA, GB>(a.G() + b.G(),
                           a.H() + b.H());
}

template <typename T, typename GA, typename GB>
Pairing<GA, GB> operator* (const T& a, const Pairing<GA, GB>& b) {
    return Pairing<GA, GB>(a * b.G(),
                           a * b.H());
}

template <typename GA, typename GB>
Pairing<GA, GB> fastAddSpecial(const Pairing<GA, GB>& a,
                               const Pairing<GA, GB>& b) {
    return Pairing<GA, GB>(fastAddSpecial(a.G(), b.G()),
                           fastAddSpecial(a.H(), b.H()));
}

template <typename GA, typename GB>
SparseVector<Pairing<GA, GB>>& batchSpecial(SparseVector<Pairing<GA, GB>>& vec)
{
    std::vector<GA> G_vec(vec.size());
    std::vector<GB> H_vec(vec.size());

    for (std::size_t i = 0; i < vec.size(); ++i) {
        G_vec[i] = vec.getElement(i).G();
        H_vec[i] = vec.getElement(i).H();
    }

    batchSpecial(G_vec);
    batchSpecial(H_vec);

    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec.setElement(i,
                       Pairing<GA, GB>(G_vec[i], H_vec[i]));
    }

    return vec;
}

template <mp_size_t N, typename GA, typename GB>
Pairing<GA, GB> wnafExp(const BigInt<N>& scalar,
                        const Pairing<GA, GB>& base)
{
    return Pairing<GA, GB>(wnafExp(scalar, base.G()),
                           wnafExp(scalar, base.H()));
}

// standard vector, works with map-reduce or monolithic window tables
template <typename GA, typename GB, typename FR>
SparseVector<Pairing<GA, GB>> batchExp(const WindowExp<GA>& tableA,
                                       const WindowExp<GB>& tableB,
                                       const FR& coeffA,
                                       const FR& coeffB,
                                       const std::vector<FR>& vec,
                                       ProgressCallback* callback = nullptr)
{
    const std::size_t M = callback ? callback->minorSteps() : 0;
    const std::size_t N = vec.size();

    SparseVector<Pairing<GA, GB>> res(N, Pairing<GA, GB>::zero());

    std::size_t i = 0, index = 0;

    // full blocks
    for (std::size_t j = 0; j < M; ++j) {
        for (std::size_t k = 0; k < N / M; ++k) {
            if (! vec[i].isZero()) {
                res.setIndexElement(
                    index++,
                    i,
                    Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                    tableB.exp(coeffB * vec[i])));
            }

            ++i;
        }

        callback->minor();
    }

    // remaining steps smaller than one block
    while (i < N) {
        if (! vec[i].isZero()) {
            res.setIndexElement(
                index++,
                i,
                Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                tableB.exp(coeffB * vec[i])));
        }

        ++i;
    }

    res.resize(index);

#ifdef USE_ADD_SPECIAL
    batchSpecial(res);
#endif

    return res;
}

// block partitioned vector, works with map-reduce or monolithic window table
template <typename GA, typename GB, typename FR>
SparseVector<Pairing<GA, GB>> batchExp(const WindowExp<GA>& tableA,
                                       const WindowExp<GB>& tableB,
                                       const FR& coeffA,
                                       const FR& coeffB,
                                       const BlockVector<FR>& vec)
{
    SparseVector<Pairing<GA, GB>> res(vec.size(), Pairing<GA, GB>::zero());

    std::size_t index = 0;

    for (std::size_t i = vec.startIndex(); i < vec.stopIndex(); ++i) {
        if (! vec[i].isZero())
            res.setIndexElement(
                index++,
                i,
                Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                tableB.exp(coeffB * vec[i])));
    }

    res.resize(index);

#ifdef USE_ADD_SPECIAL
    batchSpecial(res);
#endif

    return res;
}

// used with map-reduce
template <typename GA, typename GB, typename FR>
void batchExp(SparseVector<Pairing<GA, GB>>& res, // returned from batchExp()
              const WindowExp<GA>& tableA,
              const WindowExp<GB>& tableB,
              const FR& coeffA,
              const FR& coeffB,
              const std::vector<FR>& vec)
{
    // iterate over sparse vector directly
    for (std::size_t index = 0; index < res.size(); ++index) {
        const auto i = res.getIndex(index);
        const auto elem = res.getElement(index);

        res.setIndexElement(
            index,
            i,
            elem + Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                   tableB.exp(coeffB * vec[i])));
    }
}

// block partitioned vector, used with map-reduce
template <typename GA, typename GB, typename FR>
void batchExp(SparseVector<Pairing<GA, GB>>& res, //returned from batchExp()
              const WindowExp<GA>& tableA,
              const WindowExp<GB>& tableB,
              const FR& coeffA,
              const FR& coeffB,
              const BlockVector<FR>& vec)
{
    // iterate over sparse vector directly
    for (std::size_t index = 0; index < res.size(); ++index) {
        const auto i = res.getIndex(index);
        const auto elem = res.getElement(index);

        res.setIndexElement(
            index,
            i,
            elem + Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                   tableB.exp(coeffB * vec[i])));
    }
}

template <typename GA, typename GB, typename FR>
Pairing<GA, GB> multiExp01(const SparseVector<Pairing<GA, GB>>& base,
                           const std::vector<FR>& scalar,
                           const std::size_t minIndex,
                           const std::size_t maxIndex,
                           ProgressCallback* callback = nullptr)
{
    const auto
        ZERO = FR::zero(),
        ONE = FR::one();

    std::vector<Pairing<GA, GB>> base2;
    std::vector<FR> scalar2;

    auto accum = Pairing<GA, GB>::zero();

    for (std::size_t i = 0; i < base.size(); ++i) {
        const auto idx = base.getIndex(i);

        if (idx >= maxIndex) {
            break;

        } else if (idx >= minIndex) {
            const auto a = scalar[idx - minIndex];

            if (ZERO == a) {
                continue;

            } else if (ONE == a) {
#ifdef USE_ADD_SPECIAL
                accum = fastAddSpecial(accum, base.getElement(i));
#else
                accum = accum + base.getElement(i);
#endif

            } else {
                base2.emplace_back(base.getElement(i));
                scalar2.emplace_back(a);
            }
        }
    }

    return accum + multiExp(base2, scalar2, callback);
}

} // namespace snarklib

#endif
