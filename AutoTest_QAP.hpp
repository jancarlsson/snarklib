#ifndef _SNARKLIB_AUTOTEST_QAP_HPP_
#define _SNARKLIB_AUTOTEST_QAP_HPP_

#include <cstdint>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "qap/qap.hpp"
#include /*libsnark*/ "r1cs/r1cs.hpp"
#else
#include /*libsnark*/ "relations/arithmetic_programs/qap/qap.hpp"
#include /*libsnark*/ "relations/constraint_satisfaction_problems/r1cs/r1cs.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/AutoTest_R1CS.hpp"
#include "snarklib/QAP_query.hpp"
#include "snarklib/QAP_witness.hpp"
#include "snarklib/Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// QAP instance map
//

template <template <typename> class SYS, typename T, typename U>
class AutoTest_QAP_ABCH_instance_map : public AutoTest
{
public:
    AutoTest_QAP_ABCH_instance_map(const AutoTestR1CS<SYS, T, U>& cs,
                                   const T& point)
        : AutoTest(cs, point),
          m_constraintSystem(cs),
          m_B(point)
    {
        copyData(m_B, m_A);
    }

    AutoTest_QAP_ABCH_instance_map(const AutoTestR1CS<SYS, T, U>& cs)
        : AutoTest_QAP_ABCH_instance_map{cs, T::random()}
    {}

    void runTest() {
#ifdef USE_OLD_LIBSNARK
        const auto qapA
            = libsnark::qap_instance_map(m_constraintSystem.systemA(),
                                         m_A);
#else
        const auto qapA
            = libsnark::r1cs_to_qap_instance_map_with_evaluation(
                  m_constraintSystem.systemA(),
                  m_A);
#endif

        const QAP_SystemPoint<SYS, T> qapB(m_constraintSystem.systemB(),
                                           m_constraintSystem.numCircuitInputs(),
                                           m_B);

        if (resultsMatch(qapA, qapB)) {
            checkPass(true);

        } else {
            // try again with A and B swapped
            // libsnark swaps A and B during key generation
            // snarklib swaps A and B before key generation
            auto copyCS = m_constraintSystem.systemB();
            copyCS.swap_AB();

            const QAP_SystemPoint<SYS, T> qapB(copyCS,
                                               m_constraintSystem.numCircuitInputs(),
                                               m_B);

            checkPass(resultsMatch(qapA, qapB));
        }
    }

private:
    template <typename A>
    bool resultsMatch(const A& qapA, const QAP_SystemPoint<SYS, T>& qapB) {
        const QAP_QueryABC<SYS, T> ABCt(qapB);
        const QAP_QueryH<SYS, T> Ht(qapB);

#ifdef USE_OLD_LIBSNARK
        return
            sameData(qapA.At, ABCt.vecA()) &&
            sameData(qapA.Bt, ABCt.vecB()) &&
            sameData(qapA.Ct, ABCt.vecC()) &&
            sameData(qapA.Ht, Ht.vec()) &&
            (qapA.non_zero_At == ABCt.nonzeroA()) &&
            (qapA.non_zero_Bt == ABCt.nonzeroB()) &&
            (qapA.non_zero_Ct == ABCt.nonzeroC()) &&
            (qapA.non_zero_Ht == Ht.nonzeroCount());
#else
        // inhomogeneous terms in QAP for ABC are separate in new libsnark
        // old libsnark stored the inhomogeneous terms in the first three
        // vector elements, followed by homogeneous terms
        // snarklib follows the convention of old libsnark
        const auto
            lenA = qapA.At.size(),
            lenB = qapA.Bt.size(),
            lenC = qapA.Ct.size();

        if ((lenA + 3 != ABCt.vecA().size()) ||
            (lenB + 3 != ABCt.vecB().size()) ||
            (lenC + 3 != ABCt.vecC().size()) ||
            (lenA != lenB) ||
            (lenA != lenC)) return false;

        // compare the ABC homogeneous terms only
        for (std::size_t i = 0; i < lenA; ++i) {
            if (!sameData(qapA.At[i], ABCt.vecA()[i + 3]) ||
                !sameData(qapA.Bt[i], ABCt.vecB()[i + 3]) ||
                !sameData(qapA.Ct[i], ABCt.vecC()[i + 3])) return false;
        }

        // compare H
        return sameData(qapA.Ht, Ht.vec());
#endif
    }

    AutoTestR1CS<SYS, T, U> m_constraintSystem;
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// QAP witness map
//

template <template <typename> class SYS, typename T, typename U>
class AutoTest_QAP_Witness_map : public AutoTest
{
public:
    AutoTest_QAP_Witness_map(const AutoTestR1CS<SYS, T, U>& cs,
                             const T& d1,
                             const T& d2,
                             const T& d3)
        : AutoTest(cs, d1, d2, d3),
          m_constraintSystem(cs),
          m_d1B(d1),
          m_d2B(d2),
          m_d3B(d3)
    {
        copyData(m_d1B, m_d1A);
        copyData(m_d2B, m_d2A);
        copyData(m_d3B, m_d3A);
    }

    AutoTest_QAP_Witness_map(const AutoTestR1CS<SYS, T, U>& cs)
        : AutoTest_QAP_Witness_map{cs, T::random(), T:random(), T::random()}
    {}

    void runTest() {
#ifdef USE_OLD_LIBSNARK
        const auto HA
            = libsnark::qap_witness_map(m_constraintSystem.systemA(),
                                        m_constraintSystem.witnessA(),
                                        m_d1A,
                                        m_d2A,
                                        m_d3A);
#else
        const auto HA
            = libsnark::r1cs_to_qap_witness_map(m_constraintSystem.systemA(),
                                                m_constraintSystem.inputA(),
                                                m_constraintSystem.witnessA(),
                                                m_d1A,
                                                m_d2A,
                                                m_d3A);
#endif

        const QAP_SystemPoint<SYS, T> qap(m_constraintSystem.systemB(),
                                          m_constraintSystem.numCircuitInputs());

        const auto& HB = witnessH(qap);

#ifdef USE_OLD_LIBSNARK
        if (sameData(HA, HB)) {
#else
        if (sameData(HA.coefficients_for_H, HB)) {
#endif
            checkPass(true);

        } else {
            // try again with A and B swapped
            // libsnark swaps A and B during key generation
            // snarklib swaps A and B before key generation
            auto copyCS = m_constraintSystem.systemB();
            copyCS.swap_AB();

            const QAP_SystemPoint<SYS, T> qap(copyCS,
                                              m_constraintSystem.numCircuitInputs());

            const auto& HB = witnessH(qap);

#ifdef USE_OLD_LIBSNARK
            checkPass(sameData(HA, HB));
#else
            checkPass(sameData(HA.coefficients_for_H, HB));
#endif
        }
    }

private:
    std::vector<T> witnessH(const QAP_SystemPoint<SYS, T>& qap) const {
        const QAP_WitnessABCH<SYS, T> ABCH(qap,
                                           m_constraintSystem.witnessB(),
                                           m_d1B,
                                           m_d2B,
                                           m_d3B);

        return ABCH.vec();
    }

    AutoTestR1CS<SYS, T, U> m_constraintSystem;
    U m_d1A, m_d2A, m_d3A;
    const T m_d1B, m_d2B, m_d3B;
};

} // namespace snarklib

#endif
