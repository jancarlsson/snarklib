#ifndef _SNARKLIB_AUTOTEST_WINDOW_EXP_HPP_
#define _SNARKLIB_AUTOTEST_WINDOW_EXP_HPP_

#include <cstdint>
#include <vector>
#include "AutoTest.hpp"
#include "encoding/multiexp.hpp"
#include "WindowExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// window table size matches original
//

template <typename T, typename U>
class AutoTest_WindowExpSize : public AutoTest
{
public:
    AutoTest_WindowExpSize(const std::size_t exp_count)
        : AutoTest(exp_count),
          m_exp_count(exp_count)
    {}

    void runTest() {
        const auto a = libsnark::get_exp_window_size<U>(m_exp_count);
        const auto b = WindowExp<T>::windowBits(m_exp_count);

        checkPass(a == b);
    }

private:
    const std::size_t m_exp_count;
};

////////////////////////////////////////////////////////////////////////////////
// window table exponentiation matches original
//

template <typename T, typename F, typename U, typename G>
class AutoTest_WindowExp_exp : public AutoTest
{
public:
    AutoTest_WindowExp_exp(const std::size_t exp_count,
                           const F& value)
        : AutoTest(exp_count),
          m_exp_count(exp_count),
          m_B(value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_WindowExp_exp(const std::size_t exp_count)
        : AutoTest_WindowExp_exp{exp_count, F::random()}
    {}

    void runTest() {
        const auto a = libsnark::get_exp_window_size<U>(m_exp_count);
        const auto b = WindowExp<T>::windowBits(m_exp_count);
        if (! checkPass(a == b)) return;

        const auto A = libsnark::get_window_table(G::num_bits, U::zero(), a, U::one());
        WindowExp<T> B(b);

        const auto valueA = libsnark::windowed_exp(G::num_bits, a, A, m_A);
        const auto valueB = B.exp(m_B);

        checkPass(sameData(valueA, valueB));
    }

private:
    const std::size_t m_exp_count;
    G m_A;
    const F m_B;
};

////////////////////////////////////////////////////////////////////////////////
// window table batch exponentiation matches original
//

template <typename T, typename F, typename U, typename G>
class AutoTest_WindowExp_batchExp : public AutoTest
{
public:
    AutoTest_WindowExp_batchExp(const std::size_t exp_count,
                                const std::size_t vecSize)
        : AutoTest(exp_count),
          m_exp_count(exp_count),
          m_vecSize(vecSize),
          m_A(vecSize, G::zero())
    {
        m_B.reserve(vecSize);
        for (std::size_t i = 0; i < vecSize; ++i) {
            m_B.emplace_back(F::random());
            copyData(m_B[i], m_A[i]);
        }
    }

    void runTest() {
        const auto a = libsnark::get_exp_window_size<U>(m_exp_count);
        const auto b = WindowExp<T>::windowBits(m_exp_count);
        if (! checkPass(a == b)) return;

        const auto A = libsnark::get_window_table(G::num_bits, U::zero(), a, U::one());
        WindowExp<T> B(b);

        const auto valueA = libsnark::batch_exp(G::num_bits, a, A, m_A);
        const auto valueB = B.batchExp(m_B);

        if (checkPass(valueA.size() == valueB.size()) &&
            checkPass(valueA.size() == m_vecSize))
        {
            for (std::size_t i = 0; i < m_vecSize; ++i){
                checkPass(sameData(m_A[i], m_B[i]));
            }
        }
    }

private:
    const std::size_t m_exp_count, m_vecSize;
    std::vector<G> m_A;
    std::vector<F> m_B;
};

} // namespace snarklib

#endif
