#ifndef _SNARKLIB_QAP_WITNESS_HPP_
#define _SNARKLIB_QAP_WITNESS_HPP_

#include <cassert>
#include <cstdint>
#include <vector>
#include "HugeSystem.hpp"
#include "ProgressCallback.hpp"
#include "QAP_system.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// witness vectors A, B, C
//

template <template <typename> class SYS, typename T>
class QAP_WitnessABC
{
public:
    QAP_WitnessABC(const QAP_SystemPoint<SYS, T>& qap,
                   const R1Witness<T>& witness)
        : m_qap(qap),
          m_witness(witness),
          m_vecA(qap.degree(), T::zero()),
          m_vecB(qap.degree(), T::zero()),
          m_vecC(qap.degree(), T::zero()),
          m_uitA(m_vecA.begin()),
          m_uitB(m_vecB.begin()),
          m_uitC(m_vecC.begin()),
          m_error(false)
    {
        // input consistency (for A only)
        *m_uitA = T::one();
        for (std::size_t i = 0; i < qap.numCircuitInputs(); ++i) {
            *m_uitA += witness[i] * T(i + 2);
        }

        // A, B, C
        constraintLoop(qap.constraintSystem());

        qap.FFT()->iFFT(m_vecA);
        qap.FFT()->iFFT(m_vecB);
        qap.FFT()->iFFT(m_vecC);
    }

    void cosetFFT() {
        m_qap.FFT()->cosetFFT(m_vecA, T::params.multiplicative_generator());
        m_qap.FFT()->cosetFFT(m_vecB, T::params.multiplicative_generator());
        m_qap.FFT()->cosetFFT(m_vecC, T::params.multiplicative_generator());
    }

    const std::vector<T>& vecA() const { return m_vecA; }
    const std::vector<T>& vecB() const { return m_vecB; }
    const std::vector<T>& vecC() const { return m_vecC; }

    bool operator! () const { return m_error; }

private:
    void constraintLoop(const R1System<T>& S) {
        for (const auto& constraint : S.constraints()) {
            ++m_uitA; *m_uitA += constraint.a() * m_witness;
            ++m_uitB; *m_uitB += constraint.b() * m_witness;
            ++m_uitC; *m_uitC += constraint.c() * m_witness;
        }
    }

    void constraintLoop(const HugeSystem<T>& S) {
        m_error = S.mapLambda(
            [this] (const R1System<T>& a) -> bool {
                this->constraintLoop(a);
                return false; // do not write back to disk
            });
    }

    const QAP_SystemPoint<SYS, T>& m_qap;
    const R1Witness<T>& m_witness;
    std::vector<T> m_vecA, m_vecB, m_vecC;
    typename std::vector<T>::iterator m_uitA, m_uitB, m_uitC;
    bool m_error;
};

////////////////////////////////////////////////////////////////////////////////
// witness vector H subsumes ABC
//

template <template <typename> class SYS, typename T>
class QAP_WitnessABCH
{
public:
    QAP_WitnessABCH(const QAP_SystemPoint<SYS, T>& qap,
                    const R1Witness<T>& witness,
                    const T& random_d1,
                    const T& random_d2,
                    const T& random_d3,
                    ProgressCallback* callback = nullptr)
        : m_vec(qap.degree() + 1, T::zero())
    {
        const std::size_t N = qap.degree();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        QAP_WitnessABC<SYS, T> ABC(qap, witness);

        std::size_t i = 0;

        // for full blocks
        for (std::size_t j = 0; j < M / 2; ++j) {
            for (std::size_t k = 0; k < N / (M / 2); ++k) {
                m_vec[i] =
                    random_d2 * ABC.vecA()[i] +
                    random_d1 * ABC.vecB()[i];
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < N) {
            m_vec[i] =
                random_d2 * ABC.vecA()[i] +
                random_d1 * ABC.vecB()[i];
            ++i;
        }

        m_vec[0] -= random_d3;

        qap.FFT()->add_poly_Z(random_d1 * random_d2,
                              m_vec);

        ABC.cosetFFT();

        addTemporary(
            QAP_WitnessABCH<SYS, T>(qap, ABC),
            callback);
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    // temporary H
    QAP_WitnessABCH(const QAP_SystemPoint<SYS, T>& qap,
                    const QAP_WitnessABC<SYS, T>& ABC) // after cosetFFT()
        : m_vec(qap.degree(), T::zero())
    {
        for (std::size_t i = 0; i < m_vec.size(); ++i)
            m_vec[i] =
                ABC.vecA()[i] * ABC.vecB()[i] -
                ABC.vecC()[i];

        qap.FFT()->divide_by_Z_on_coset(m_vec);

        qap.FFT()->icosetFFT(m_vec,
                             T::params.multiplicative_generator());
    }

    // add regular and temporary H together
    void addTemporary(const QAP_WitnessABCH& tmpH,
                      ProgressCallback* callback = nullptr)
    {
        const std::size_t N = tmpH.vec().size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

#ifdef USE_ASSERT
        // make sure to add temporary H, not regular H
        assert(N < m_vec.size());
#endif

        std::size_t i = 0;

        // for full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                m_vec[i] += tmpH.vec()[i];
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < N) {
            m_vec[i] += tmpH.vec()[i];
            ++i;
        }
    }

    std::vector<T> m_vec;
};

} // namespace snarklib

#endif
