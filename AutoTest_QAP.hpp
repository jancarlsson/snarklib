#ifndef _SNARKLIB_AUTOTEST_QAP_HPP_
#define _SNARKLIB_AUTOTEST_QAP_HPP_

#include <cstdint>
#include "AutoTest.hpp"
#include "AutoTest_R1CS.hpp"
#include "QAP_query.hpp"
#include "QAP_witness.hpp"
#include "qap/qap.hpp"
#include "Rank1DSL.hpp"
#include "r1cs/r1cs.hpp"

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
        const auto qapA
            = libsnark::qap_instance_map(m_constraintSystem.systemA(),
                                         m_A);

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

        return
            sameData(qapA.At, ABCt.vecA()) &&
            sameData(qapA.Bt, ABCt.vecB()) &&
            sameData(qapA.Ct, ABCt.vecC()) &&
            sameData(qapA.Ht, Ht.vec()) &&
            (qapA.non_zero_At == ABCt.nonzeroA()) &&
            (qapA.non_zero_Bt == ABCt.nonzeroB()) &&
            (qapA.non_zero_Ct == ABCt.nonzeroC()) &&
            (qapA.non_zero_Ht == Ht.nonzeroCount());
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
        const auto HA
            = libsnark::qap_witness_map(m_constraintSystem.systemA(),
                                        m_constraintSystem.witnessA(),
                                        m_d1A,
                                        m_d2A,
                                        m_d3A);

        const QAP_SystemPoint<SYS, T> qap(m_constraintSystem.systemB(),
                                          m_constraintSystem.numCircuitInputs());

        const auto& HB = witnessH(qap);

        if (sameData(HA, HB)) {
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

            checkPass(sameData(HA, HB));
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
