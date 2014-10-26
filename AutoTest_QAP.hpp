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

        const QAP_ABCH<T> qapB(m_constraintSystem.systemB(),
                               m_constraintSystem.numberInputs(),
                               m_B);

        checkPass(sameData(qapA.At, qapB.A_query()));
        checkPass(sameData(qapA.Bt, qapB.B_query()));
        checkPass(sameData(qapA.Ct, qapB.C_query()));
        checkPass(sameData(qapA.Ht, qapB.H_query()));

        checkPass(qapA.non_zero_At == qapB.nonzeroAt());
        checkPass(qapA.non_zero_Bt == qapB.nonzeroBt());
        checkPass(qapA.non_zero_Ct == qapB.nonzeroCt());
        checkPass(qapA.non_zero_Ht == qapB.nonzeroHt());
    }

private:
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

        const QAP_Witness<T> qapB(m_constraintSystem.systemB(),
                                  m_constraintSystem.numberInputs(),
                                  m_constraintSystem.witnessB(),
                                  m_d1B,
                                  m_d2B,
                                  m_d3B);

        const auto& HB = qapB.H();

        checkPass(sameData(HA, HB));
    }

private:
    const AutoTestR1CS<T, U> m_constraintSystem;
    U m_d1A, m_d2A, m_d3A;
    const T m_d1B, m_d2B, m_d3B;
};

} // namespace snarklib

#endif
