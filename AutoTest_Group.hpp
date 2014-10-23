#ifndef _SNARKLIB_AUTOTEST_GROUP_HPP_
#define _SNARKLIB_AUTOTEST_GROUP_HPP_

#include <gmp.h>
#include <string>
#include "AutoTest.hpp"
#include "algebra/fields/bigint.hpp"
#include "BigInt.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// zero remains zero
//

template <typename T>
class AutoTest_GroupZeroStaysZero : public AutoTest
{
public:
    AutoTest_GroupZeroStaysZero()
        : AutoTest()
    {}

    void runTest() {
        const auto ZERO = T::zero();

        const auto
            a = ZERO + ZERO,
            b = ZERO - ZERO,
            c = BigInt<1>().randomize() * ZERO,
            d = BigInt<2>().randomize() * ZERO,
            e = BigInt<3>().randomize() * ZERO;

        const auto f = fastAddSpecial(a, b);

        checkPass(ZERO == a && a == ZERO && a.isZero());
        checkPass(ZERO == b && b == ZERO && b.isZero());
        checkPass(ZERO == c && c == ZERO && c.isZero());
        checkPass(ZERO == d && d == ZERO && d.isZero());
        checkPass(ZERO == e && e == ZERO && e.isZero());
        checkPass(ZERO == f && f == ZERO && f.isZero());
    }
};

////////////////////////////////////////////////////////////////////////////////
// additive identity
//

template <typename T>
class AutoTest_GroupZeroIdentity : public AutoTest
{
public:
    AutoTest_GroupZeroIdentity(const T& value)
        : AutoTest(value),
          m_B(value)
    {}

    AutoTest_GroupZeroIdentity()
        : AutoTest_GroupZeroIdentity{T::random()}
    {}

    void runTest() {
        const auto
            a = m_B + T::zero(),
            b = T::zero() + m_B,
            c = m_B - T::zero();

        const auto
            d = fastAddSpecial(a, T::zero()),
            e = fastAddSpecial(T::zero(), b);

        checkPass(a == b && a == m_B && b == m_B && c == m_B && m_B == c);
        checkPass(d == e && d == a && e == b);
    }

private:
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// multiplicative identity
//

template <typename T>
class AutoTest_GroupOneIdentity : public AutoTest
{
    typedef typename T::ScalarField Fr;

public:
    AutoTest_GroupOneIdentity(const T& value)
        : AutoTest(value),
          m_B(value)
    {}

    AutoTest_GroupOneIdentity()
        : AutoTest_GroupOneIdentity{T::random()}
    {}

    void runTest() {
        const auto
            a = Fr::random(),
            b = Fr::random();

        checkPass((a + b) * T::one() == a * T::one() + b * T::one());

        const auto
            c = a - Fr::one(),
            d = a + Fr::one();

        checkPass((c + Fr::one()) * T::one() == (d - Fr::one()) * T::one());
        checkPass(c * T::one() == a * T::one() - T::one());
        checkPass(d * T::one() == a * T::one() + T::one());

        auto
            specA = a * T::one(),
            specB = b * T::one(),
            specAB = (a + b) * T::one();

        specA.toSpecial();
        specB.toSpecial();
        specAB.toSpecial();

        checkPass(specAB == fastAddSpecial(specA, specB));
    }

private:
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// addition matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_GroupAdd : public AutoTest
{
public:
    AutoTest_GroupAdd(const std::string& lhs, const std::string& rhs)
        : AutoTest(lhs, rhs),
          m_lhsA(to_bigint<N>(lhs) * U::one()),
          m_rhsA(to_bigint<N>(rhs) * U::one()),
          m_lhsB(BigInt<N>(lhs) * T::one()),
          m_rhsB(BigInt<N>(rhs) * T::one())
    {}

    void runTest() {
        const auto a = m_lhsA + m_rhsA;
        const auto b = m_lhsB + m_rhsB;
        checkPass(sameData(a, b));

        auto
            a1 = m_lhsA,
            a2 = m_rhsA;
        a1.to_special();
        a2.to_special();

        auto
            b1 = m_lhsB,
            b2 = m_rhsB;
        b1.toSpecial();
        b2.toSpecial();

        checkPass(
            sameData(
                a1.fast_add_special(a2),
                fastAddSpecial(b1, b2)));
    }

private:
    const U m_lhsA, m_rhsA;
    const T m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// subtraction matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_GroupSub : public AutoTest
{
public:
    AutoTest_GroupSub(const std::string& lhs, const std::string& rhs)
        : AutoTest(lhs, rhs),
          m_lhsA(to_bigint<N>(lhs) * U::one()),
          m_rhsA(to_bigint<N>(rhs) * U::one()),
          m_lhsB(BigInt<N>(lhs) * T::one()),
          m_rhsB(BigInt<N>(rhs) * T::one())
    {}

    void runTest() {
        const auto a = m_lhsA - m_rhsA;
        const auto b = m_lhsB - m_rhsB;
        checkPass(sameData(a, b));

        auto
            a1 = m_lhsA,
            a2 = -m_rhsA;
        a1.to_special();
        a2.to_special();

        auto
            b1 = m_lhsB,
            b2 = -m_rhsB;
        b1.toSpecial();
        b2.toSpecial();

        checkPass(
            sameData(
                a1.fast_add_special(a2),
                fastAddSpecial(b1, b2)));
    }

private:
    const U m_lhsA, m_rhsA;
    const T m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// multiplication matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_GroupMul : public AutoTest
{
public:
    AutoTest_GroupMul(const std::string& pow, const std::string& base)
        : AutoTest(pow, base),
          m_powerA(pow.c_str()),
          m_powerB(pow),
          m_baseA(to_bigint<N>(base) * U::one()),
          m_baseB(BigInt<N>(base) * T::one())
    {}

    void runTest() {
        const auto a = m_powerA * m_baseA;
        const auto b = m_powerB * m_baseB;

        checkPass(sameData(a, b));
    }

private:
    const libsnark::bigint<N> m_powerA;
    const BigInt<N> m_powerB;
    const U m_baseA;
    const T m_baseB;
};

////////////////////////////////////////////////////////////////////////////////
// doubling matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_GroupDbl : public AutoTest
{
public:
    AutoTest_GroupDbl(const std::string& value)
        : AutoTest(value),
          m_A(to_bigint<N>(value) * U::one()),
          m_B(BigInt<N>(value) * T::one())
    {}

    void runTest() {
        const auto a = m_A.dbl();
        const auto b = m_B.dbl();
        
        checkPass(sameData(a, b));
    }

private:
    const U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// special and well formed matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_GroupSpecialWellFormed : public AutoTest
{
public:
    AutoTest_GroupSpecialWellFormed(const std::string& value)
        : AutoTest(value),
          m_A(to_bigint<N>(value) * U::one()),
          m_B(BigInt<N>(value) * T::one())
    {}

    void runTest() {
        auto a = m_A;
        auto b = m_B;

        for (std::size_t i = 0; i < 1000; ++i) {
            checkPass(a.is_well_formed() == b.wellFormed());
            checkPass(a.is_special() == b.isSpecial());

            auto c = a;
            auto d = b;
            c.to_special();
            d.toSpecial();

            checkPass(c.is_well_formed() == d.wellFormed());
            checkPass(c.is_special() == d.isSpecial());

            a = a + U::one();
            b = b + T::one();
        }
    }

private:
    const U m_A;
    const T m_B;
};

} // namespace snarklib

#endif
