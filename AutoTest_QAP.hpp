#ifndef _SNARKLIB_AUTOTEST_QAP_HPP_
#define _SNARKLIB_AUTOTEST_QAP_HPP_

#include <cstdint>
#include "AutoTest.hpp"
#include "AutoTest_R1CS.hpp"
#include "QAP.hpp"
#include "qap/qap.hpp"
#include "Rank1DSL.hpp"
#include "r1cs/r1cs.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// QAP instance map
//

template <typename T, typename U>
class AutoTest_QAP_ABCH_instance_map : public AutoTest
{
public:
    AutoTest_QAP_ABCH_instance_map(const AutoTestR1CS<T, U>& cs,
                                   const T& point)
        : AutoTest(cs, point),
          m_constraintSystem(cs),
          m_B(point)
    {
        copyData(m_B, m_A);
    }

    AutoTest_QAP_ABCH_instance_map(const AutoTestR1CS<T, U>& cs)
        : AutoTest_QAP_ABCH_instance_map{cs, T::random()}
    {}

    void runTest() {
        const auto qapA
            = libsnark::qap_instance_map(m_constraintSystem.systemA(),
                                         m_A);

        const QAP_SystemPoint<T> qapB(m_constraintSystem.systemB(),
                                      m_constraintSystem.numberInputs(),
                                      m_B);

        if (resultsMatch(qapA, qapB)) {
            checkPass(true);

        } else {
            // try again with A and B swapped
            // libsnark swaps A and B during key generation
            // snarklib swaps A and B before key generation
            auto copyCS = m_constraintSystem.systemB();
            copyCS.swap_AB();

            const QAP_SystemPoint<T> qapB(copyCS,
                                          m_constraintSystem.numberInputs(),
                                          m_B);

            checkPass(resultsMatch(qapA, qapB));
        }
    }

private:
    template <typename A>
    bool resultsMatch(const A& qapA, const QAP_SystemPoint<T>& qapB) {
        const QAP_QueryA<T> At(qapB);
        const QAP_QueryB<T> Bt(qapB);
        const QAP_QueryC<T> Ct(qapB);
        const QAP_QueryH<T> Ht(qapB);

        return
            sameData(qapA.At, At.vec()) &&
            sameData(qapA.Bt, Bt.vec()) &&
            sameData(qapA.Ct, Ct.vec()) &&
            sameData(qapA.Ht, Ht.vec()) &&
            (qapA.non_zero_At == At.nonzeroCount()) &&
            (qapA.non_zero_Bt == Bt.nonzeroCount()) &&
            (qapA.non_zero_Ct == Ct.nonzeroCount()) &&
            (qapA.non_zero_Ht == Ht.nonzeroCount());
    }

    const AutoTestR1CS<T, U> m_constraintSystem;
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// QAP witness map
//

template <typename T, typename U>
class AutoTest_QAP_Witness_map : public AutoTest
{
public:
    AutoTest_QAP_Witness_map(const AutoTestR1CS<T, U>& cs,
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

    AutoTest_QAP_Witness_map(const AutoTestR1CS<T, U>& cs)
        : AutoTest_QAP_Witness_map{cs, T::random(), T:random(), T::random()}
    {}

    void runTest() {
        const auto HA
            = libsnark::qap_witness_map(m_constraintSystem.systemA(),
                                        m_constraintSystem.witnessA(),
                                        m_d1A,
                                        m_d2A,
                                        m_d3A);

        const QAP_SystemPoint<T> qap(m_constraintSystem.systemB(),
                                     m_constraintSystem.numberInputs());

        const auto& HB = witnessH(qap);

        if (sameData(HA, HB)) {
            checkPass(true);

        } else {
            // try again with A and B swapped
            // libsnark swaps A and B during key generation
            // snarklib swaps A and B before key generation
            auto copyCS = m_constraintSystem.systemB();
            copyCS.swap_AB();

            const QAP_SystemPoint<T> qap(copyCS,
                                         m_constraintSystem.numberInputs());

            const auto& HB = witnessH(qap);

            checkPass(sameData(HA, HB));
        }
    }

private:
    std::vector<T> witnessH(const QAP_SystemPoint<T>& qap) const {
        QAP_WitnessA<T> aA(qap, m_constraintSystem.witnessB());
        QAP_WitnessB<T> aB(qap, m_constraintSystem.witnessB());
        QAP_WitnessC<T> aC(qap, m_constraintSystem.witnessB());
        QAP_WitnessH<T> aH(qap, aA, aB, m_d1B, m_d2B, m_d3B);

        aA.cosetFFT();
        aB.cosetFFT();
        aC.cosetFFT();

        aH.addTemporary(QAP_WitnessH<T>(qap, aA, aB, aC));

        return aH.vec();
    }

    const AutoTestR1CS<T, U> m_constraintSystem;
    U m_d1A, m_d2A, m_d3A;
    const T m_d1B, m_d2B, m_d3B;
};

} // namespace snarklib

#endif
