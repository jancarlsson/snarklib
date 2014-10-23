#ifndef _SNARKLIB_PAIRING_HPP_
#define _SNARKLIB_PAIRING_HPP_

#include <cstdint>
#include <gmp.h>
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "Group.hpp"
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

    static Pairing zero() {
        return Pairing<GA, GB>(GA::zero(), GB::zero());
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
    std::vector<GA> G_vec;
    std::vector<GB> H_vec;
    G_vec.reserve(vec.size());
    H_vec.reserve(vec.size());

    for (std::size_t i = 0; i < vec.size(); ++i) {
        G_vec.push_back(vec.getElement(i).G());
        H_vec.push_back(vec.getElement(i).H());
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

template <typename GA, typename GB, typename FR>
SparseVector<Pairing<GA, GB>> batchExp(const WindowExp<GA>& tableA,
                                       const WindowExp<GB>& tableB,
                                       const FR& coeffA,
                                       const FR& coeffB,
                                       const std::vector<FR>& vec)
{
    SparseVector<Pairing<GA, GB>> res(vec.size());

    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (! vec[i].isZero())
            res.pushBack(i,
                         Pairing<GA, GB>(tableA.exp(coeffA * vec[i]),
                                         tableB.exp(coeffB * vec[i])));
    }

#ifdef USE_ADD_SPECIAL
    batchSpecial(res);
#endif

    return res;
}

template <typename GA, typename GB, typename FR>
Pairing<GA, GB> multiExp01(const SparseVector<Pairing<GA, GB>>& base,
                           const std::vector<FR>& scalar,
                           const std::size_t minIndex,
                           const std::size_t maxIndex)
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

    return accum + multiExp(base2, scalar2);
}

} // namespace snarklib

#endif
