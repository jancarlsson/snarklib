#ifndef _SNARKLIB_AUTOTEST_MULTIEXP_HPP_
#define _SNARKLIB_AUTOTEST_MULTIEXP_HPP_

#include <gmp.h>
#include <string>
#include <vector>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "algebra/fields/bigint.hpp"
#include /*libsnark*/ "common/wnaf.hpp"
#include /*libsnark*/ "encoding/multiexp.hpp"
#else
#include /*libsnark*/ "algebra/fields/bigint.hpp"
#include /*libsnark*/ "algebra/scalar_multiplication/wnaf.hpp"
#include /*libsnark*/ "algebra/scalar_multiplication/multiexp.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/BigInt.hpp"
#include "snarklib/MultiExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// wNAF exponentiation matches original
//

template <mp_size_t N, typename T, typename U>
class AutoTest_MultiExp_wnafExp : public AutoTest
{
public:
    AutoTest_MultiExp_wnafExp(const std::string& scalar,
                              const std::string& base)
        : AutoTest(scalar, base),
          m_scalarA(scalar.c_str()),
          m_scalarB(scalar),
          m_baseA(to_bigint<N>(base) * U::one()),
          m_baseB(BigInt<N>(base) * T::one())
    {}

    void runTest() {
        const auto a = opt_window_wnaf_exp(
#ifdef USE_OLD_LIBSNARK
                                           U::zero(),
#endif
                                           m_baseA,
                                           m_scalarA,
                                           m_scalarA.num_bits());

        const auto b = wnafExp(m_scalarB, m_baseB);

        checkPass(sameData(a, b));
    }

private:
    const libsnark::bigint<N> m_scalarA;
    const BigInt<N> m_scalarB;
    const U m_baseA;
    const T m_baseB;
};

////////////////////////////////////////////////////////////////////////////////
// sum of multiple exponentiation matches original
//

template <mp_size_t N, typename T, typename F, typename U, typename G>
class AutoTest_MultiExp_multiExp : public AutoTest
{
public:
    AutoTest_MultiExp_multiExp(const std::size_t numTerms)
        : AutoTest(numTerms),
          m_numTerms(numTerms)
    {
        m_baseB.reserve(numTerms);
        m_scalarB.reserve(numTerms);
        m_baseA.reserve(numTerms);
        m_scalarA.reserve(numTerms);

        for (std::size_t i = 0; i < numTerms; ++i) {
            const auto
                randomBase = uniformBase10(0, 1000000),
                randomScalar = uniformBase10(0, 1000000);

            m_baseB.emplace_back(T(BigInt<N>(randomBase) * T::one()));
            m_scalarB.emplace_back(F(randomScalar));

            m_baseA.emplace_back(U(to_bigint<N>(randomBase) * U::one()));
            m_scalarA.emplace_back(G(randomScalar.c_str()));
        }
    }

    void runTest() {
        const auto a = libsnark::multi_exp<U, G>(
#ifdef USE_OLD_LIBSNARK
                                                 U::zero(),
#endif
                                                 m_baseA.begin(),
                                                 m_baseA.end(),
                                                 m_scalarA.begin(),
                                                 m_scalarA.end(),
                                                 1,
                                                 true);

        const auto b = multiExp(m_baseB, m_scalarB);

        checkPass(sameData(a, b));
    }

private:
    const std::size_t m_numTerms;
    std::vector<T> m_baseB;
    std::vector<F> m_scalarB;
    std::vector<U> m_baseA;
    std::vector<G> m_scalarA;
};

////////////////////////////////////////////////////////////////////////////////
// sum of multiple exponentiation with zeros and ones matches original
//

template <mp_size_t N, typename T, typename F, typename U, typename G>
class AutoTest_MultiExp_multiExp01 : public AutoTest
{
public:
    AutoTest_MultiExp_multiExp01(const std::size_t numTerms)
        : AutoTest(numTerms),
          m_numTerms(numTerms)
    {
        m_baseB.reserve(numTerms);
        m_scalarB.reserve(numTerms);
        m_baseA.reserve(numTerms);
        m_scalarA.reserve(numTerms);

        for (std::size_t i = 0; i < numTerms; ++i) {
            const auto
                randomBase = uniformBase10(0, 1000000),
                randomScalar = sparseUniformBase10(0, 1000000);

            m_baseB.emplace_back(T(BigInt<N>(randomBase) * T::one()));
            m_scalarB.emplace_back(F(randomScalar));

            m_baseA.emplace_back(U(to_bigint<N>(randomBase) * U::one()));
            m_scalarA.emplace_back(G(randomScalar.c_str()));
        }
    }

    void runTest() {
        const auto a = libsnark::multi_exp<U, G>(
#ifdef USE_OLD_LIBSNARK
                                                 U::zero(),
#endif
                                                 m_baseA.begin(),
                                                 m_baseA.end(),
                                                 m_scalarA.begin(),
                                                 m_scalarA.end(),
                                                 1,
                                                 true);

        const auto b = multiExp(m_baseB, m_scalarB);

        checkPass(sameData(a, b));
    }

private:
    const std::size_t m_numTerms;
    std::vector<T> m_baseB;
    std::vector<F> m_scalarB;
    std::vector<U> m_baseA;
    std::vector<G> m_scalarA;
};

} // namespace snarklib

#endif
