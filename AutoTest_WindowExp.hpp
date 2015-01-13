#ifndef _SNARKLIB_AUTOTEST_WINDOW_EXP_HPP_
#define _SNARKLIB_AUTOTEST_WINDOW_EXP_HPP_

#include <cstdint>
#include <vector>
#include "AutoTest.hpp"
#include "AuxSTL.hpp"
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
        : AutoTest(exp_count, value),
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
        const WindowExp<T> B(m_exp_count);

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
        : AutoTest(exp_count, vecSize),
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
        const WindowExp<T> B(m_exp_count);

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

////////////////////////////////////////////////////////////////////////////////
// compare map-reduce with monolithic window table exponentiation
//

template <typename T, typename F>
class AutoTest_WindowExp_expMapReduce : public AutoTest
{
public:
    AutoTest_WindowExp_expMapReduce(const std::size_t exp_count,
                                    const F& value)
        : AutoTest(exp_count, value),
          m_exp_count(exp_count),
          m_value(value)
    {}

    AutoTest_WindowExp_expMapReduce(const std::size_t exp_count)
        : AutoTest_WindowExp_expMapReduce{exp_count, F::random()}
    {}

    void runTest() {
        const WindowExp<T> A(m_exp_count);
        const auto result_A = A.exp(m_value);

        const auto space = WindowExp<T>::space(m_exp_count);

        // try all possible block partitionings
        for (std::size_t numBlocks = 1; numBlocks <= space.globalID()[0]; ++numBlocks) {
            auto idx = space;
            idx.blockPartition(std::array<std::size_t, 1>{ numBlocks });

            auto result_B = T::zero();

            // mapping
            for (std::size_t block = 0; block < numBlocks; ++block) {
                const WindowExp<T> B(idx, block);

                // reducing
                result_B = result_B + B.exp(m_value);
            }

            checkPass(result_A == result_B);
        }
    }

private:
    const std::size_t m_exp_count;
    const F m_value;
};

////////////////////////////////////////////////////////////////////////////////
// compare block partitioned with standard vector batch exponentiation
//

template <typename T, typename F>
class AutoTest_WindowExp_batchExpMapReduce1 : public AutoTest
{
public:
    AutoTest_WindowExp_batchExpMapReduce1(const std::size_t exp_count,
                                          const std::size_t vecSize)
        : AutoTest(exp_count, vecSize),
          m_exp_count(exp_count),
          m_vecSize(vecSize)
    {
        m_vec.reserve(vecSize);
        for (std::size_t i = 0; i < vecSize; ++i)
            m_vec.emplace_back(F::random());
    }

    void runTest() {
        const WindowExp<T> A(m_exp_count);
        const auto result_A = A.batchExp(m_vec);

        const auto space = BlockVector<F>::space(m_vec);

        // try all possible block partitionings
        for (std::size_t numBlocks = 1; numBlocks <= space.globalID()[0]; ++numBlocks) {
            auto idx = space;
            idx.blockPartition(std::array<std::size_t, 1>{ numBlocks });

            std::vector<T> result_B(m_vecSize);

            for (std::size_t block = 0; block < numBlocks; ++block) {
                BlockVector<F> partvec(idx, block, m_vec);
                A.batchExp(partvec).emplace(result_B);
            }

            checkPass(result_A == result_B);
        }
    }

private:
    const std::size_t m_exp_count, m_vecSize;
    std::vector<F> m_vec;
};

////////////////////////////////////////////////////////////////////////////////
// map-reduce window tables and block partitioned vector batch exponentiation
//

template <typename T, typename F>
class AutoTest_WindowExp_batchExpMapReduce2 : public AutoTest
{
public:
    AutoTest_WindowExp_batchExpMapReduce2(const std::size_t exp_count,
                                          const std::size_t vecSize)
        : AutoTest(exp_count, vecSize),
          m_exp_count(exp_count),
          m_vecSize(vecSize)
    {
        m_vec.reserve(vecSize);
        for (std::size_t i = 0; i < vecSize; ++i)
            m_vec.emplace_back(F::random());
    }

    void runTest() {
        const WindowExp<T> A(m_exp_count);
        const auto result_A = A.batchExp(m_vec);

        const auto winSpace = WindowExp<T>::space(m_exp_count);
        const auto vecSpace = BlockVector<F>::space(m_vec);

        // just try three partitionings of window table
        for (const auto numWinBlks : std::array<std::size_t, 3>{ 1, 2, winSpace.globalID()[0] }) {
            auto winIdx = winSpace;
            winIdx.blockPartition(std::array<std::size_t, 1>{ numWinBlks });

            // try all possible block partitionings of vector
            for (std::size_t numVecBlks = 1; numVecBlks <= vecSpace.globalID()[0]; ++numVecBlks) {
                auto vecIdx = vecSpace;
                vecIdx.blockPartition(std::array<std::size_t, 1>{ numVecBlks });

                std::vector<T> result_B(m_vecSize);

                // (outer loop) iterate over windows
                for (std::size_t winblock = 0; winblock < numWinBlks; ++winblock) {
                    // partial window table is expensive
                    // ***must be in outer loop***
                    const WindowExp<T> partwin(winIdx, winblock);

                    // (inner loop) iterate over vector blocks
                    for (std::size_t vecblock = 0; vecblock < numVecBlks; ++vecblock) {
                        // read in
                        BlockVector<T> result(vecIdx, vecblock, result_B);

                        // accumulate from partial window table
                        const BlockVector<F> partvec(vecIdx, vecblock, m_vec);
                        result += partwin.batchExp(partvec);

                        // write back
                        result.emplace(result_B);
                    }
                }

                checkPass(result_A == result_B);
            }
        }
    }

private:
    const std::size_t m_exp_count, m_vecSize;
    std::vector<F> m_vec;
};

} // namespace snarklib

#endif
