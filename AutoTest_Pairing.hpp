#ifndef _SNARKLIB_AUTOTEST_PAIRING_HPP_
#define _SNARKLIB_AUTOTEST_PAIRING_HPP_

#include <algorithm>
#include <gmp.h>
#include <random>
#include <string>
#include <vector>
#include "algebra/fields/bigint.hpp"
#include "AutoTest.hpp"
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "encoding/knowledge_commitment.hpp"
#include "encoding/multiexp.hpp"
#include "Pairing.hpp"
#include "WindowExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// paired groups addition matches original
//

template <mp_size_t N, typename TG, typename TH, typename UG, typename UH>
class AutoTest_PairingAdd : public AutoTest
{
public:
    AutoTest_PairingAdd(const std::string& a,
                        const std::string& b,
                        const std::string& c,
                        const std::string& d)
        : AutoTest(a, b, c, d),
          m_A1({to_bigint<N>(a) * UG::one(), to_bigint<N>(b) * UH::one()}),
          m_A2({to_bigint<N>(c) * UG::one(), to_bigint<N>(d) * UH::one()}),
          m_B1(BigInt<N>(a) * TG::one(), BigInt<N>(b) * TH::one()),
          m_B2(BigInt<N>(c) * TG::one(), BigInt<N>(d) * TH::one())
    {}

    void runTest() {
        const auto a = m_A1 + m_A2;
        const auto b = m_B1 + m_B2;

        checkPass(sameData(a, b));
    }

private:
    const libsnark::knowledge_commitment<UG, UH> m_A1, m_A2;
    const Pairing<TG, TH> m_B1, m_B2;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups multiplication matches original
//

template <mp_size_t N,
          typename TG, typename TH, typename TF,
          typename UG, typename UH, typename UF>
class AutoTest_PairingMul : public AutoTest
{
public:
    AutoTest_PairingMul(const TF& a,
                        const std::string& b,
                        const std::string& c)
        : AutoTest(a, b, c),
          m_fB(a),
          m_gA({to_bigint<N>(b) * UG::one(), to_bigint<N>(c) * UH::one()}),
          m_gB(BigInt<N>(b) * TG::one(), BigInt<N>(c) * TH::one())
    {
        copyData(m_fB, m_fA);
    }

    AutoTest_PairingMul(const std::string& b,
                        const std::string& c)
        : AutoTest_PairingMul{TF::random(), b, c}
    {}

    void runTest() {
        const auto a = m_fA * m_gA;
        const auto b = m_fB * m_gB;

        checkPass(sameData(a, b));
    }

private:
    UF m_fA;
    const TF m_fB;
    const libsnark::knowledge_commitment<UG, UH> m_gA;
    const Pairing<TG, TH> m_gB;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups fast addition with zeros and ones matches original
//

template <mp_size_t N, typename TG, typename TH, typename UG, typename UH>
class AutoTest_PairingFastAddSpecial : public AutoTest
{
public:
    AutoTest_PairingFastAddSpecial(const std::string& a,
                                   const std::string& b,
                                   const std::string& c,
                                   const std::string& d)
        : AutoTest(a, b, c, d),
          m_A1({to_bigint<N>(a) * UG::one(), to_bigint<N>(b) * UH::one()}),
          m_A2({to_bigint<N>(c) * UG::one(), to_bigint<N>(d) * UH::one()}),
          m_B1(BigInt<N>(a) * TG::one(), BigInt<N>(b) * TH::one()),
          m_B2(BigInt<N>(c) * TG::one(), BigInt<N>(d) * TH::one())
    {}

    void runTest() {
        m_A1.g = m_A1.g.fast_add_special(m_A2.g);
        m_A1.h = m_A1.h.fast_add_special(m_A2.h);

        const auto b = fastAddSpecial(m_B1, m_B2);

        checkPass(sameData(m_A1, b));
    }

private:
    libsnark::knowledge_commitment<UG, UH> m_A1;
    const libsnark::knowledge_commitment<UG, UH> m_A2;
    const Pairing<TG, TH> m_B1, m_B2;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups sparse vector to special representation matches original
//

template <mp_size_t N, typename TG, typename TH, typename UG, typename UH>
class AutoTest_PairingBatchSpecial : public AutoTest
{
public:
    AutoTest_PairingBatchSpecial(const std::size_t vecSize)
        : AutoTest(vecSize),
          m_vecSize(vecSize),
          m_B(vecSize, Pairing<TG, TH>::zero())
    {
        m_A.values.reserve(vecSize);
        m_A.indices.reserve(vecSize);

        std::random_device rd;
        std::size_t idx = 0;

        for (std::size_t i = 0; i < vecSize; ++i) {
            const std::string
                a = randomBase10(rd, N),
                b = randomBase10(rd, N);

            m_A.indices.push_back(idx);
            m_A.values.emplace_back(
                libsnark::knowledge_commitment<UG, UH>(
                    to_bigint<N>(a) * UG::one(), to_bigint<N>(b) * UH::one()));
            m_A.is_sparse = true;
            m_A.original_size = 0;

            m_B.setIndexElement(
                i,
                idx,
                Pairing<TG, TH>(BigInt<N>(a) * TG::one(), BigInt<N>(b) * TH::one()));

            idx += 1 + rd() % 10;
        }
    }

    void runTest() {
        libsnark::kc_batch_to_special(m_A.values);
        batchSpecial(m_B);

        if (checkPass(m_A.values.size() == m_vecSize) &&
            checkPass(m_A.indices.size() == m_vecSize) &&
            checkPass(m_B.size() == m_vecSize))
        {
            for (std::size_t i = 0; i < m_vecSize; ++i) {
                checkPass(m_A.indices[i] == m_B.getIndex(i));
                checkPass(sameData(m_A.values[i], m_B.getElement(i)));
            }
        }
    }

private:
    const std::size_t m_vecSize;
    libsnark::knowledge_commitment_vector<UG, UH> m_A;
    SparseVector<Pairing<TG, TH>> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups wnafExp matches original
//

template <mp_size_t N, typename TG, typename TH, typename UG, typename UH>
class AutoTest_Pairing_wnafExp : public AutoTest
{
public:
    AutoTest_Pairing_wnafExp(const std::string& a,
                             const std::string& b,
                             const std::string& c)
        : AutoTest(a, b, c),
          m_scalarA(a.c_str()),
          m_scalarB(a),
          m_baseA({to_bigint<N>(b) * UG::one(), to_bigint<N>(c) * UH::one()}),
          m_baseB(BigInt<N>(b) * TG::one(), BigInt<N>(c) * TH::one())
    {}

    void runTest() {
        const auto a = libsnark::opt_window_wnaf_exp(
            libsnark::knowledge_commitment<UG, UH>(UG::zero(), UH::zero()),
            m_baseA,
            m_scalarA,
            m_scalarA.num_bits());

        const auto b = wnafExp(m_scalarB, m_baseB);

        checkPass(sameData(a, b));
    }

private:
    const libsnark::bigint<N> m_scalarA;
    const BigInt<N> m_scalarB;
    const libsnark::knowledge_commitment<UG, UH> m_baseA;
    const Pairing<TG, TH> m_baseB;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups batchExp matches original
//

template <mp_size_t N,
          typename TG, typename TH, typename TF,
          typename UG, typename UH, typename UF>
class AutoTest_Pairing_batchExp : public AutoTest
{
public:
    AutoTest_Pairing_batchExp(const std::size_t exp_count,
                              const TF& coeffA,
                              const TF& coeffB,
                              const std::size_t vecSize)
        : AutoTest(exp_count, coeffA, coeffB, vecSize),
          m_exp_count(exp_count),
          m_coeffA_B(coeffA),
          m_coeffB_B(coeffB),
          m_vecA(vecSize, UF::zero())
    {
        copyData(m_coeffA_B, m_coeffA_A);
        copyData(m_coeffB_B, m_coeffB_A);

        m_vecB.reserve(vecSize);
        for (std::size_t i = 0; i < vecSize; ++i) {
            // note: zero element makes kc_batch_exp() crash
            m_vecB.emplace_back(TF(uniformBase10(1, 1000000)));
            copyData(m_vecB[i], m_vecA[i]);
        }
    }

    AutoTest_Pairing_batchExp(const std::size_t exp_count,
                              const std::size_t vecSize)
        : AutoTest_Pairing_batchExp{exp_count, TF::random(), TF::random(), vecSize}
    {}

    void runTest() {
        // exp_count to window table sizes
        const auto
            expSize_UG_A = libsnark::get_exp_window_size<UG>(m_exp_count),
            expSize_UH_A = libsnark::get_exp_window_size<UH>(m_exp_count);

        const auto
            expSize_TG_B = WindowExp<TG>::windowBits(m_exp_count),
            expSize_TH_B = WindowExp<TH>::windowBits(m_exp_count);

        if (! checkPass(expSize_UG_A == expSize_TG_B) ||
            ! checkPass(expSize_UH_A == expSize_TH_B)) {
            return;
        }

        // create window tables
        const auto winTable_UG_A = libsnark::get_window_table(UF::num_bits,
                                                              UG::zero(),
                                                              expSize_UG_A,
                                                              UG::one());

        const auto winTable_UH_A = libsnark::get_window_table(UF::num_bits,
                                                              UH::zero(),
                                                              expSize_UH_A,
                                                              UH::one());

        const WindowExp<TG> winTable_TG_B(m_exp_count);
        const WindowExp<TH> winTable_TH_B(m_exp_count);

        // batch exponentiation
        const auto sparseVec_A = kc_batch_exp(UF::num_bits,
                                              expSize_UG_A,
                                              expSize_UH_A,
                                              winTable_UG_A,
                                              winTable_UH_A,
                                              m_coeffA_A,
                                              m_coeffB_A,
                                              m_vecA,
                                              true,
                                              1);

        const auto sparseVec_B = batchExp(winTable_TG_B,
                                          winTable_TH_B,
                                          m_coeffA_B,
                                          m_coeffB_B,
                                          m_vecB);

        // compare sparse vector output
        if (checkPass(sparseVec_A.values.size() == sparseVec_B.size()) &&
            checkPass(sparseVec_A.indices.size() == sparseVec_B.size()))
        {
            for (std::size_t i = 0; i < sparseVec_B.size(); ++i) {
                checkPass(sparseVec_A.indices[i] == sparseVec_B.getIndex(i));
                checkPass(sameData(sparseVec_A.values[i], sparseVec_B.getElement(i)));
            }
        }
    }

private:
    const std::size_t m_exp_count;
    UF m_coeffA_A, m_coeffB_A;
    const TF m_coeffA_B, m_coeffB_B;
    std::vector<UF> m_vecA;
    std::vector<TF> m_vecB;
};

////////////////////////////////////////////////////////////////////////////////
// paired groups multiExp01 matches original
//

template <mp_size_t N,
          typename TG, typename TH, typename TF,
          typename UG, typename UH, typename UF>
class AutoTest_Pairing_multiExp01 : public AutoTest
{
public:
    AutoTest_Pairing_multiExp01(const std::size_t vecSize)
        : AutoTest(vecSize),
          m_baseB(vecSize, Pairing<TG, TH>::zero())
    {
        m_baseA.values.reserve(vecSize);
        m_baseA.indices.reserve(vecSize);

        std::random_device rd;
        std::size_t idx = 0;

        for (std::size_t i = 0; i < vecSize; ++i) {
            const std::string
                a = randomBase10(rd, N),
                b = randomBase10(rd, N);

            m_baseA.indices.push_back(idx);
            m_baseA.values.emplace_back(
                libsnark::knowledge_commitment<UG, UH>(
                    to_bigint<N>(a) * UG::one(), to_bigint<N>(b) * UH::one()));
            m_baseA.is_sparse = true;
            m_baseA.original_size = 0;

            m_baseB.setIndexElement(
                i,
                idx,
                Pairing<TG, TH>(BigInt<N>(a) * TG::one(), BigInt<N>(b) * TH::one()));

            idx += 1 + rd() % 10;
        }

        // arbitrary choice of minIndex and maxIndex inside [0, idx]
        m_minIndex = idx / 10;
        m_maxIndex = idx / 2;

        const std::size_t scalarSize = m_maxIndex - m_minIndex + 1;
        m_scalarA.reserve(scalarSize);
        m_scalarB.reserve(scalarSize);
        for (std::size_t i = 0; i < scalarSize; ++i) {
            const std::string a = sparseUniformBase10(0, 1000000);
            m_scalarB.emplace_back(TF(a));
            m_scalarA.emplace_back(UF(a.c_str()));
            copyData(m_scalarB[i], m_scalarA[i]);
        }
    }

    void runTest() {
        const auto a = libsnark::kc_multi_exp_with_fast_add_special<UG, UH, UF>(
            libsnark::knowledge_commitment<UG, UH>(UG::zero(), UH::zero()),
            m_baseA,
            m_minIndex,
            m_maxIndex,
            m_scalarA.begin(),
            m_scalarA.end(),
            1,
            true);

        const auto b = multiExp01(m_baseB, m_scalarB, m_minIndex, m_maxIndex);

        checkPass(sameData(a, b));
    }

private:
    libsnark::knowledge_commitment_vector<UG, UH> m_baseA;
    SparseVector<Pairing<TG, TH>> m_baseB;
    std::vector<UF> m_scalarA;
    std::vector<TF> m_scalarB;
    std::size_t m_minIndex, m_maxIndex;
};

////////////////////////////////////////////////////////////////////////////////
// compare map-reduce with monolithic batchExp
//

template <typename TG, typename TH, typename TF>
class AutoTest_Pairing_batchExpMapReduce1 : public AutoTest
{
public:
    AutoTest_Pairing_batchExpMapReduce1(const std::size_t exp_count,
                                        const TF& coeffA,
                                        const TF& coeffB,
                                        const std::size_t vecSize)
        : AutoTest(exp_count, coeffA, coeffB, vecSize),
          m_exp_count(exp_count),
          m_coeffA(coeffA),
          m_coeffB(coeffB),
          m_vec(vecSize, TF::zero())
    {
        for (std::size_t i = 0; i < vecSize; ++i)
            m_vec[i] = TF(uniformBase10(1, 1000000));
    }

    AutoTest_Pairing_batchExpMapReduce1(const std::size_t exp_count,
                                        const std::size_t vecSize)
        : AutoTest_Pairing_batchExpMapReduce1{exp_count, TF::random(), TF::random(), vecSize}
    {}

    void runTest() {
        const WindowExp<TG> tableG(m_exp_count);
        const WindowExp<TH> tableH(m_exp_count);

        const auto result_A = batchExp(tableG, tableH, m_coeffA, m_coeffB, m_vec);

        const auto space = WindowExp<TG>::space(m_exp_count);

        // try all possible block partitionings
        for (std::size_t numBlocks = 1; numBlocks <= space.globalID()[0]; ++numBlocks) {
            auto idx = space;
            idx.blockPartition(std::array<std::size_t, 1>{ numBlocks });

            // initial block
            const WindowExp<TG> startG(idx, 0);
            const WindowExp<TH> startH(idx, 0);
            auto result_B = batchExp(startG, startH, m_coeffA, m_coeffB, m_vec);

            // mapping subsequent blocks
            for (std::size_t block = 1; block < numBlocks; ++block) {
                const WindowExp<TG> partG(idx, block);
                const WindowExp<TH> partH(idx, block);

                // reducing subsequent blocks
                batchExp(result_B, partG, partH, m_coeffA, m_coeffB, m_vec);
            }

            checkPass(result_A == result_B);
        }
    }

private:
    const std::size_t m_exp_count;
    const TF m_coeffA, m_coeffB;
    std::vector<TF> m_vec;
};

////////////////////////////////////////////////////////////////////////////////
// map-reduce window tables and block partitioned vector batchExp
//

template <typename TG, typename TH, typename TF>
class AutoTest_Pairing_batchExpMapReduce2 : public AutoTest
{
public:
    AutoTest_Pairing_batchExpMapReduce2(const std::size_t exp_count,
                                        const TF& coeffA,
                                        const TF& coeffB,
                                        const std::size_t vecSize)
        : AutoTest(exp_count, coeffA, coeffB, vecSize),
          m_exp_count(exp_count),
          m_coeffA(coeffA),
          m_coeffB(coeffB),
          m_vec(vecSize, TF::zero())
    {
        for (std::size_t i = 0; i < vecSize; ++i)
            m_vec[i] = TF(uniformBase10(1, 1000000));
    }

    AutoTest_Pairing_batchExpMapReduce2(const std::size_t exp_count,
                                        const std::size_t vecSize)
        : AutoTest_Pairing_batchExpMapReduce2{exp_count, TF::random(), TF::random(), vecSize}
    {}

    void runTest() {
        const WindowExp<TG> tableG(m_exp_count);
        const WindowExp<TH> tableH(m_exp_count);

        const auto result_A = batchExp(tableG, tableH, m_coeffA, m_coeffB, m_vec);

        const auto winSpace = WindowExp<TG>::space(m_exp_count);
        const auto vecSpace = BlockVector<TF>::space(m_vec);

        // just try three partitionings of window table
        for (const auto numWinBlks : std::array<std::size_t, 3>{ 1, 2, winSpace.globalID()[0] }) {
            auto winIdx = winSpace;
            winIdx.blockPartition(std::array<std::size_t, 1>{ numWinBlks });

            // try all possible block partitionings of vector
            for (std::size_t numVecBlks = 1; numVecBlks <= vecSpace.globalID()[0]; ++numVecBlks) {
                auto vecIdx = vecSpace;
                vecIdx.blockPartition(std::array<std::size_t, 1>{ numVecBlks });

                std::vector<SparseVector<Pairing<TG, TH>>> result(numVecBlks);

                // initial block
                const WindowExp<TG> startG(winIdx, 0);
                const WindowExp<TH> startH(winIdx, 0);
                for (std::size_t vecblock = 0; vecblock < numVecBlks; ++vecblock) {
                    BlockVector<TF> partvec(vecIdx, vecblock, m_vec);
                    result[vecblock] = batchExp(startG, startH, m_coeffA, m_coeffB, partvec);
                }

                // subsequent blocks
                for (std::size_t winblock = 1; winblock < numWinBlks; ++winblock) {
                    const WindowExp<TG> partG(winIdx, winblock);
                    const WindowExp<TH> partH(winIdx, winblock);
                    for (std::size_t vecblock = 0; vecblock < numVecBlks; ++vecblock) {
                        BlockVector<TF> partvec(vecIdx, vecblock, m_vec);
                        batchExp(result[vecblock], partG, partH, m_coeffA, m_coeffB, partvec);
                    }
                }

                auto& result_B = result[0];
                for (std::size_t i = 1; i < numVecBlks; ++i)
                    result_B.concat(result[i]);

                checkPass(result_A == result_B);
            }
        }
    }

private:
    const std::size_t m_exp_count;
    const TF m_coeffA, m_coeffB;
    std::vector<TF> m_vec;
};

} // namespace snarklib

#endif
