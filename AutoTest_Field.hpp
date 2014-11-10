#ifndef _SNARKLIB_AUTOTEST_FIELD_HPP_
#define _SNARKLIB_AUTOTEST_FIELD_HPP_

#include <gmp.h>
#include <string>
#include "AutoTest.hpp"
#include "algebra/fields/bigint.hpp"
#include "BigInt.hpp"
#include "FpX.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// zero remains zero
//

template <typename T>
class AutoTest_FieldZeroStaysZero : public AutoTest
{
public:
    AutoTest_FieldZeroStaysZero()
        : AutoTest()
    {}

    void runTest() {
        const auto ZERO = T::zero();

        const auto
            a = ZERO + ZERO,
            b = ZERO - ZERO,
            c = ZERO * ZERO,
            d = ZERO ^ 1,
            e = ZERO ^ 2,
            f = ZERO ^ 3;

        checkPass(ZERO == a && a == ZERO && a.isZero());
        checkPass(ZERO == b && b == ZERO && b.isZero());
        checkPass(ZERO == c && c == ZERO && c.isZero());
        checkPass(ZERO == d && d == ZERO && d.isZero());
        checkPass(ZERO == e && e == ZERO && e.isZero());
        checkPass(ZERO == f && f == ZERO && f.isZero());
    }
};

////////////////////////////////////////////////////////////////////////////////
// one remains one
//

template <typename T>
class AutoTest_FieldOneStaysOne : public AutoTest
{
public:
    AutoTest_FieldOneStaysOne()
        : AutoTest()
    {}

    void runTest() {
        const auto ONE = T::one();

        const auto
            a = ONE * ONE,
            b = ONE ^ 0,
            c = ONE ^ 1,
            d = ONE ^ 2,
            e = ONE ^ 3;

        checkPass(ONE == a && a == ONE);
        checkPass(ONE == b && b == ONE);
        checkPass(ONE == c && c == ONE);
        checkPass(ONE == d && d == ONE);
        checkPass(ONE == e && e == ONE);
    }
};

////////////////////////////////////////////////////////////////////////////////
// additive and multiplicative identities
//

template <typename T>
class AutoTest_FieldZeroAndOneIdentities : public AutoTest
{
public:
    AutoTest_FieldZeroAndOneIdentities(const T& value)
        : AutoTest(value),
          m_B(value)
    {}

    AutoTest_FieldZeroAndOneIdentities()
        : AutoTest_FieldZeroAndOneIdentities{T::random()}
    {}

    void runTest() {
        auto v = m_B;

        v = v + T::zero();
        v = v * T::one();
        v = T::zero() + v;
        v = T::one() * v;

        checkPass(m_B == v);
    }

private:
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// addition matches original
//

template <typename T, typename U>
class AutoTest_FieldAdd : public AutoTest
{
public:
    AutoTest_FieldAdd(const T& lhs, const T& rhs)
        : AutoTest(lhs, rhs),
          m_lhsB(lhs),
          m_rhsB(rhs)
    {
        copyData(m_lhsB, m_lhsA);
        copyData(m_rhsB, m_rhsA);
    }

    AutoTest_FieldAdd()
        : AutoTest_FieldAdd{T::random(), T::random()}
    {}

    void runTest() {
        const auto a = m_lhsA + m_rhsA;
        const auto b = m_lhsB + m_rhsB;

        checkPass(sameData(a, b));
    }

private:
    U m_lhsA, m_rhsA;
    const T m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// subtraction matches original
//

template <typename T, typename U>
class AutoTest_FieldSub : public AutoTest
{
public:
    AutoTest_FieldSub(const T& lhs, const T& rhs)
        : AutoTest(lhs, rhs),
          m_lhsB(lhs),
          m_rhsB(rhs)
    {
        copyData(m_lhsB, m_lhsA);
        copyData(m_rhsB, m_rhsA);
    }

    AutoTest_FieldSub()
        : AutoTest_FieldSub{T::random(), T::random()}
    {}

    void runTest() {
        const auto a = m_lhsA - m_rhsA;
        const auto b = m_lhsB - m_rhsB;

        checkPass(sameData(a, b));
    }

private:
    U m_lhsA, m_rhsA;
    const T m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// multiplication matches original
//

template <typename T, typename U>
class AutoTest_FieldMul : public AutoTest
{
public:
    AutoTest_FieldMul(const T& lhs, const T& rhs)
        : AutoTest(lhs, rhs),
          m_lhsB(lhs),
          m_rhsB(rhs)
    {
        copyData(m_lhsB, m_lhsA);
        copyData(m_rhsB, m_rhsA);
    }

    AutoTest_FieldMul()
        : AutoTest_FieldMul{T::random(), T::random()}
    {}

    void runTest() {
        const auto a = m_lhsA * m_rhsA;
        const auto b = m_lhsB * m_rhsB;

        checkPass(sameData(a, b));
    }

private:
    U m_lhsA, m_rhsA;
    const T m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// exponentiation with big integer matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_FieldExp : public AutoTest
{
public:
    AutoTest_FieldExp(const T& base, const BigInt<N>& pow)
        : AutoTest(base, pow),
          m_B(base),
          m_powerB(pow)
    {
        copyData(m_B, m_A);
        copyData(m_powerB, m_powerA);
    }

    AutoTest_FieldExp()
        : AutoTest_FieldExp{T::random(), BigInt<N>::random()}
    {}

    void runTest() {
        const auto a = m_A ^ m_powerA;
        const auto b = m_B ^ m_powerB;

        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
    libsnark::bigint<N> m_powerA;
    const BigInt<N> m_powerB;
};

////////////////////////////////////////////////////////////////////////////////
// squaring matches original
//

template <typename T, typename U>
class AutoTest_FieldSquared : public AutoTest
{
public:
    AutoTest_FieldSquared(const T& value)
        : AutoTest(value),
          m_B(value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_FieldSquared()
        : AutoTest_FieldSquared{T::random()}
    {}

    void runTest() {
        const auto a = m_A.squared();
        const auto b = squared(m_B);
        
        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// square root matches original
//

template <typename T, typename U>
class AutoTest_FieldSqrt : public AutoTest
{
public:
    AutoTest_FieldSqrt(const T& value)
        : AutoTest(value),
          m_B(value * value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_FieldSqrt()
        : AutoTest_FieldSqrt{T::random()}
    {}

    void runTest() {
        const auto a = m_A.sqrt();
        const auto b = sqrt(m_B);

        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// inversion matches original
//

template <typename T, typename U>
class AutoTest_FieldInverse : public AutoTest
{
public:
    AutoTest_FieldInverse(const T& value)
        : AutoTest(value),
          m_B(value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_FieldInverse()
        : AutoTest_FieldInverse{T::random()}
    {}

    void runTest() {
        const auto a = m_A.inverse();
        const auto b = inverse(m_B);

        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// Frobenius map for: Fp2, Fp3, Fp23, Fp32, Fp232 - matches original
// (all field types except Fp)
//

template <typename T, typename U>
class AutoTest_FieldFrobeniusMap : public AutoTest
{
public:
    AutoTest_FieldFrobeniusMap(const T& value, const unsigned long pow)
        : AutoTest(value, pow),
          m_B(value),
          m_power(pow)
    {
        copyData(m_B, m_A);
    }

    AutoTest_FieldFrobeniusMap(const unsigned long pow)
        : AutoTest_FieldFrobeniusMap{T::random(), pow}
    {}

    void runTest() {
        const auto a = m_A.Frobenius_map(m_power);
        const auto b = Frobenius_map(m_B, m_power);

        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
    const unsigned long m_power;
};

////////////////////////////////////////////////////////////////////////////////
// cyclotomic exponentiation matches original
// (defined for Fp32 and Fp232, the Fqk/GT fields for BN128 and Edwards)
//

template <mp_size_t N, typename T, typename U>
class AutoTest_FieldCyclotomicExp : public AutoTest
{
public:
    AutoTest_FieldCyclotomicExp(const T& value, const BigInt<N>& pow)
        : AutoTest(value),
          m_B(value),
          m_powerB(pow)
    {
        copyData(m_B, m_A);
        copyData(m_powerB, m_powerA);
    }

    AutoTest_FieldCyclotomicExp()
        : AutoTest_FieldCyclotomicExp{T::random(), BigInt<N>::random()}
    {}

    void runTest() {
        const auto a = m_A.cyclotomic_exp(m_powerA);
        const auto b = cyclotomic_exp(m_B, m_powerB);

        checkPass(sameData(a, b));
    }

private:
    U m_A;
    const T m_B;
    libsnark::bigint<N> m_powerA;
    const BigInt<N> m_powerB;
};

////////////////////////////////////////////////////////////////////////////////
// multiplication by 024 matches original
// used by BN128 pairing only
// for arguments: (Fp232, Fp2, Fp2, Fp2)
//           i.e. (Fqk, Fqe, Fqe, Fqe)
//

template <typename TFQK, typename TFQE, typename UFQK, typename UFQE>
class AutoTest_FieldMulBy024 : public AutoTest
{
public:
    AutoTest_FieldMulBy024(const TFQK& x,
                           const TFQE& ell_0,
                           const TFQE& ell_VW,
                           const TFQE& ell_VV)
        : AutoTest(x, ell_0, ell_VW, ell_VV),
          m_xB(x),
          m_ell_0B(ell_0),
          m_ell_VWB(ell_VW),
          m_ell_VVB(ell_VV)
    {
        copyData(m_xB, m_xA);
        copyData(m_ell_0B, m_ell_0A);
        copyData(m_ell_VWB, m_ell_VWA);
        copyData(m_ell_VVB, m_ell_VVA);
    }

    AutoTest_FieldMulBy024()
        : AutoTest_FieldMulBy024{TFQK::random(),
                                 TFQE::random(),
                                 TFQE::random(),
                                 TFQE::random()}
    {}

    void runTest() {
        const auto a = m_xA.mul_by_024(m_ell_0A, m_ell_VWA, m_ell_VVA);
        const auto b = mul_by_024(m_xB, m_ell_0B, m_ell_VWB, m_ell_VVB);

        checkPass(sameData(a, b));
    }

private:
    UFQK m_xA;
    UFQE m_ell_0A, m_ell_VWA, m_ell_VVA;
    const TFQK m_xB;
    const TFQE m_ell_0B, m_ell_VWB, m_ell_VVB;
};

} // namespace snarklib

#endif
