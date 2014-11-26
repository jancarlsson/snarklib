#ifndef _SNARKLIB_GROUP_HPP_
#define _SNARKLIB_GROUP_HPP_

#include <algorithm>
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <tuple>
#include <vector>
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "FpModel.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Symmetric pairing group for elliptic curves
//

template <typename BASE, typename SCALAR, typename CURVE>
class Group
{
public:
    typedef BASE BaseField;
    typedef SCALAR ScalarField;

    // symmetric pairing group parameters
    class Params
    {
    public:
        // G_zero
        const Group& G_zero() const {
            return m_G_zero;
        }
        void G_zero(const BASE& x, const BASE& y) {
            m_G_zero = Group(x, y);
        }
        void G_zero(const BASE& x, const BASE& y, const BASE& z) {
            m_G_zero = Group(x, y, z);
        }

        // G_one
        const Group& G_one() const {
            return m_G_one;
        }
        void G_one(const BASE& x, const BASE& y) {
            m_G_one = Group(x, y);
        }
        void G_one(const BASE& x, const BASE& y, const BASE& z) {
            m_G_one = Group(x, y, z);
        }

        // wnaf_window_table
        const std::vector<std::size_t>& wnaf_window_table() const {
            return m_wnaf_window_table;
        }
        void wnaf_window_table_clear() {
            m_wnaf_window_table.clear();
        }
        void wnaf_window_table(const std::size_t a) {
            m_wnaf_window_table.push_back(a);
        }

        // fixed_base_exp_window_table
        const std::vector<std::size_t>& fixed_base_exp_window_table() const {
            return m_fixed_base_exp_window_table;
        }
        void fixed_base_exp_window_table_clear() {
            m_fixed_base_exp_window_table.clear();
        }
        void fixed_base_exp_window_table(const std::size_t a) {
            m_fixed_base_exp_window_table.push_back(a);
        }

    private:
        static Group m_G_zero;
        static Group m_G_one;

        std::vector<std::size_t> m_wnaf_window_table;
        std::vector<std::size_t> m_fixed_base_exp_window_table;
    };

    static Params params;

    // default is zero
    Group()
        : Group{zero().m_X, zero().m_Y, zero().m_Z}
    {}

    // inverted coordinates
    Group(const BASE& x, const BASE& y)
        : Group{y, x, x*y}
    {}

    // Jacobian coordinates
    Group(const BASE& x, const BASE& y, const BASE& z)
        : m_X(x), m_Y(y), m_Z(z)
    {}

    // break encapsulation (used for pairing)
    const BASE& x() const { return m_X; }
    const BASE& y() const { return m_Y; }
    const BASE& z() const { return m_Z; }
    void x(const BASE& a) { m_X = a; }
    void y(const BASE& a) { m_Y = a; }
    void z(const BASE& a) { m_Z = a; }

    void affineCoordinates() {
        std::tie(m_X, m_Y, m_Z)
            = CURVE::affineCoordinates(m_X, m_Y, m_Z);
    }

    void toSpecial() {
        std::tie(m_X, m_Y, m_Z)
            = CURVE::toSpecial(m_X, m_Y, m_Z);
    }

    bool isSpecial() const {
        return isZero() || m_Z == BASE::one();
    }

    bool isZero() const {
        return CURVE::isZero(m_X, m_Y, m_Z);
    }

    bool operator== (const Group& other) const {
        return CURVE::equalOp(
            m_X, m_Y, m_Z,
            other.m_X, other.m_Y, other.m_Z);
    }

    bool operator!= (const Group& other) const {
        return ! (*this == other);
    }

    Group operator- () const {
        return CURVE::negateOp(m_X, m_Y, m_Z, *this);
    }

    Group dbl() const {
        return CURVE::dbl(m_X, m_Y, m_Z, *this);
    }

    bool wellFormed() const {
        return CURVE::wellFormed(m_X, m_Y, m_Z);
    }

    static const Group& zero() {
        return params.G_zero();
    }

    static const Group& one() {
        return params.G_one();
    }

    static Group random() {
        return SCALAR::BaseType::random() * one();
    }

    static std::size_t sizeInBits() {
        return BASE::sizeInBits() + 1;
    }

    static constexpr
    const BigInt<BASE::BaseType::numberLimbs()>& baseModulus() {
        return BASE::BaseType::modulus();
    }

    static constexpr
    const BigInt<SCALAR::BaseType::numberLimbs()>& scalarModulus() {
        return SCALAR::BaseType::modulus();
    }

    void marshal_out(std::ostream& os) const {
        x().marshal_out(os);
        y().marshal_out(os);
        z().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_X.marshal_in(is) &&
            m_Y.marshal_in(is) &&
            m_Z.marshal_in(is);
    }

private:
    BASE m_X, m_Y, m_Z;
};

////////////////////////////////////////////////////////////////////////////////
// static data members
//

template <typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE>
Group<BASE, SCALAR, CURVE>::Params::m_G_zero;

template <typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE>
Group<BASE, SCALAR, CURVE>::Params::m_G_one;

template <typename BASE, typename SCALAR, typename CURVE>
typename Group<BASE, SCALAR, CURVE>::Params
Group<BASE, SCALAR, CURVE>::params;

////////////////////////////////////////////////////////////////////////////////
// Operator functions
//

template <typename BASE, typename SCALAR, typename CURVE>
void marshal_out(std::ostream& os,
                 const std::vector<Group<BASE, SCALAR, CURVE>>& a) {
    // size
    os << a.size() << std::endl;

    // group vector
    for (const auto& g : a) {
        g.marshal_out(os);
    }
}

template <typename BASE, typename SCALAR, typename CURVE>
bool marshal_in(std::istream& is,
                std::vector<Group<BASE, SCALAR, CURVE>>& a) {
    // size
    std::size_t numberElems;
    is >> numberElems;
    if (!is) return false;

    // group vector
    a.clear();
    a.reserve(numberElems);
    for (std::size_t i = 0; i < numberElems; ++i) {
        Group<BASE, SCALAR, CURVE> g;
        if (!g.marshal_in(is)) return false;
        a.emplace_back(g);
    }

    return true; // ok
}

template <typename BASE, typename SCALAR, typename CURVE>
std::ostream& operator<< (std::ostream& out,
                          const Group<BASE, SCALAR, CURVE>& a) {
    auto copy(a);
    copy.affineCoordinates();

    CURVE::outputPrefix(out, a);

    return out << copy.x() << " "
               << (copy.y()[0].asUnsignedLong() & 1);
}

template <typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> operator+ (const Group<BASE, SCALAR, CURVE>& a,
                                      const Group<BASE, SCALAR, CURVE>& b) {
    return CURVE::addOp(
        a.x(), a.y(), a.z(),
        b.x(), b.y(), b.z(),
        a);
}
    
template <typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> operator- (const Group<BASE, SCALAR, CURVE>& a,
                                      const Group<BASE, SCALAR, CURVE>& b) {
    return a + (-b);
}

template <typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> fastAddSpecial(const Group<BASE, SCALAR, CURVE>& a,
                                          const Group<BASE, SCALAR, CURVE>& b) {
    return CURVE::fastAddSpecial(
        a.x(), a.y(), a.z(),
        b.x(), b.y(), b.z(),
        a);
}
    
template <mp_size_t N,
          typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> operator* (const BigInt<N>& exponent,
                                      const Group<BASE, SCALAR, CURVE>& base) {
    return power(exponent, base); // group version: base follows power
                                  // this uses dbl() and operator+
}
    
template <mp_size_t N, const BigInt<N>& MODULUS,
          typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> operator* (const FpModel<N, MODULUS>& exponent,
                                      const Group<BASE, SCALAR, CURVE>& base) {
    return exponent.asBigInt() * base;
}

template <mp_size_t N, const BigInt<N>& MODULUS,
          typename BASE, typename SCALAR, typename CURVE>
Group<BASE, SCALAR, CURVE> operator* (const Field<FpModel<N, MODULUS>>& exponent,
                                      const Group<BASE, SCALAR, CURVE>& base) {
    return exponent[0] * base;
}

// batch conversion to special (batch_invert() makes it faster)
template <typename BASE, typename SCALAR, typename CURVE>
std::vector<Group<BASE, SCALAR, CURVE>>&
batchSpecial(std::vector<Group<BASE, SCALAR, CURVE>>& vec) {
    return CURVE::batchSpecial(vec);
}

} // namespace snarklib

#endif
