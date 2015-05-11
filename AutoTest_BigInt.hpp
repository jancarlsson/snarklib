#ifndef _SNARKLIB_AUTOTEST_BIGINT_HPP_
#define _SNARKLIB_AUTOTEST_BIGINT_HPP_

#include <cstdint>
#include <gmp.h>
#include <sstream>
#include <string>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "algebra/fields/bigint.hpp"
#include /*libsnark*/ "common/wnaf.hpp"
#else
#include /*libsnark*/ "algebra/fields/bigint.hpp"
#include /*libsnark*/ "algebra/scalar_multiplication/wnaf.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/BigInt.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// default constructor creates zero number
//

template <mp_size_t N>
class AutoTest_BigIntDefaultConstructorZero : public AutoTest
{
public:
    AutoTest_BigIntDefaultConstructorZero()
        : AutoTest(N)
    {}

    void runTest() {
        checkPass(m_A.is_zero() && 0 == m_A.as_ulong());
        checkPass(m_B.isZero() && 0 == m_B.asUnsignedLong());
        checkPass(sameData(m_A, m_B));
    }

private:
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// unsigned long constructor matches original
//

template <mp_size_t N>
class AutoTest_BigIntUnsignedLongConstructor : public AutoTest
{
public:
    AutoTest_BigIntUnsignedLongConstructor(const unsigned long start,
                                           const unsigned long step,
                                           const std::size_t count)
        : AutoTest(start, step, count),
          m_start(start),
          m_step(step),
          m_count(count)
    {}

    void runTest() {
        unsigned long v = m_start;

        for (std::size_t i = 0; i < m_count; ++i) {
            const libsnark::bigint<N> A(v);
            const BigInt<N> B(v);

            checkPass(sameData(A, B));
            checkPass(v == A.as_ulong());
            checkPass(v == B.asUnsignedLong());

            v += m_step;
        }
    }

private:
    const unsigned long m_start, m_step;
    const std::size_t m_count;
};

////////////////////////////////////////////////////////////////////////////////
// string constructor matches original
//

template <mp_size_t N>
class AutoTest_BigIntStringConstructor : public AutoTest
{
public:
    AutoTest_BigIntStringConstructor(const std::string& base10)
        : AutoTest(base10),
          m_A(base10.c_str()),
          m_B(base10)
    {}

    void runTest() {
        checkPass(sameData(m_A, m_B));
    }

private:
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// equality matches original
//

template <mp_size_t N>
class AutoTest_BigIntEquality : public AutoTest
{
public:
    AutoTest_BigIntEquality(const std::string& base10_LHS,
                            const std::string& base10_RHS)
        : AutoTest(base10_LHS, base10_RHS),
          m_lhsA(base10_LHS.c_str()),
          m_rhsA(base10_RHS.c_str()),
          m_lhsB(base10_LHS),
          m_rhsB(base10_RHS)
    {}

    void runTest() {
        checkPass((m_lhsA == m_rhsA) == (m_lhsB == m_rhsB));
        checkPass((m_lhsA != m_rhsA) == (m_lhsB != m_rhsB));
    }

private:
    const libsnark::bigint<N> m_lhsA, m_rhsA;
    const BigInt<N> m_lhsB, m_rhsB;
};

////////////////////////////////////////////////////////////////////////////////
// stream output matches original
//

template <mp_size_t N>
class AutoTest_BigIntStreamOutput : public AutoTest
{
public:
    AutoTest_BigIntStreamOutput(const std::string& base10)
        : AutoTest(base10),
          m_A(base10.c_str()),
          m_B(base10)
    {}

    void runTest() {
        std::stringstream ssA, ssB;
        ssA << m_A;
        ssB << m_B;

        checkPass(ssA.str() == ssB.str());
    }

private:
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// number bits matches original
//

template <mp_size_t N>
class AutoTest_BigIntNumBits : public AutoTest
{
public:
    AutoTest_BigIntNumBits(const std::string& base10)
        : AutoTest(base10),
          m_A(base10.c_str()),
          m_B(base10)
    {}

    void runTest() {
        checkPass(m_A.num_bits() == m_B.numBits());
    }

private:
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// test bits match original
//

template <mp_size_t N>
class AutoTest_BigIntTestBits : public AutoTest
{
public:
    AutoTest_BigIntTestBits(const std::string& base10)
        : AutoTest(base10),
          m_A(base10.c_str()),
          m_B(base10)
    {}

    void runTest() {
        const std::size_t
            A_bits = m_A.max_bits(),
            B_bits = m_B.maxBits();

        if (checkPass(A_bits == B_bits)) {
            for (std::size_t i = 0; i < A_bits; ++i)
                checkPass(m_A.test_bit(i) == m_B.testBit(i));
        }
    }

private:
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// weighted non-adjacent form matches original
//

template <mp_size_t N>
class AutoTest_BigIntFindwNAF : public AutoTest
{
public:
    AutoTest_BigIntFindwNAF(const std::size_t w,
                            const std::string& base10)
        : AutoTest(w, base10),
          m_w(w),
          m_A(base10.c_str()),
          m_B(base10)
    {}

    void runTest() {
#ifdef USE_OLD_LIBSNARK
        const auto a = libsnark::find_wNAF(m_w, m_A);
#else
        const auto a = libsnark::find_wnaf(m_w, m_A);
#endif
        const auto b = find_wNAF(m_w, m_B);

        if (checkPass(a.size() == b.size())) {
            for (std::size_t i = 0; i < a.size(); ++i) {
                checkPass(a[i] == b[i]);
            }
        }
    }

private:
    const std::size_t m_w;
    const libsnark::bigint<N> m_A;
    const BigInt<N> m_B;
};

} // namespace snarklib

#endif
