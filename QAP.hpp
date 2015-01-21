#ifndef _SNARKLIB_QAP_HPP_
#define _SNARKLIB_QAP_HPP_

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>
#include "LagrangeFFT.hpp"
#include "LagrangeFFTX.hpp"
#include "ProgressCallback.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// base class
//

template <typename T>
class QAP
{
public:
    virtual ~QAP() = default;

    std::size_t degree() const { return m_degree; }
    std::size_t numVariables() const { return m_numVariables; }
    std::size_t numCircuitInputs() const { return m_numCircuitInputs; }

    const LagrangeFFT<T>& FFT() const { return m_FFT; }

protected:
    QAP(const R1System<T>& constraintSystem,
        const std::size_t numCircuitInputs)
        : m_degree(LagrangeFFT<T>::getDegree(constraintSystem.constraints().size() + 1)),
          m_numVariables(constraintSystem.maxIndex()),
          m_numCircuitInputs(numCircuitInputs),
          m_FFT(m_degree)
    {}

private:
    std::size_t m_degree, m_numVariables, m_numCircuitInputs;
    LagrangeFFT<T> m_FFT;
};

////////////////////////////////////////////////////////////////////////////////
// constraint system evaluated at point
//

template <typename T>
class QAP_SystemPoint : public QAP<T>
{
public:
    // keypair generation: ABCH, IC coefficients, and K
    QAP_SystemPoint(const R1System<T>& constraintSystem,
                    const std::size_t numCircuitInputs,
                    const T& point)
        : QAP<T>(constraintSystem, numCircuitInputs),
          m_constraintSystem(constraintSystem),
          m_point(point),
          m_compute_Z(QAP<T>::FFT()->compute_Z(point)),
          m_lagrange_coeffs(QAP<T>::FFT()->lagrange_coeffs(point))
    {}

    // proof generation
    QAP_SystemPoint(const R1System<T>& constraintSystem,
                    const std::size_t numCircuitInputs)
        : QAP<T>(constraintSystem, numCircuitInputs),
          m_constraintSystem(constraintSystem),
          m_point(T::zero()),
          m_compute_Z(T::zero()),
          m_lagrange_coeffs()
    {}

    std::size_t degree() const { return QAP<T>::degree(); }
    std::size_t numVariables() const { return QAP<T>::numVariables(); }
    std::size_t numCircuitInputs() const { return QAP<T>::numCircuitInputs(); }

    const R1System<T>& constraintSystem() const { return m_constraintSystem; }
    const T& point() const { return m_point; }
    const T& compute_Z() const { return m_compute_Z; }
    const std::vector<T>& lagrange_coeffs() const { return m_lagrange_coeffs; }

private:
    const R1System<T>& m_constraintSystem;
    const T m_point, m_compute_Z;
    const std::vector<T> m_lagrange_coeffs;
};

////////////////////////////////////////////////////////////////////////////////
// query vectors A, B, C
//

template <typename T, char R1C, std::size_t Z_INDEX>
class QAP_QueryABC
{
public:
    QAP_QueryABC(const QAP_SystemPoint<T>& qap)
        : m_nonzeroCount(0),
          m_vec(3 + qap.numVariables() + 1, T::zero())
    {
        m_vec[Z_INDEX] = qap.compute_Z();
        auto uit = qap.lagrange_coeffs().begin();

        // input consistency
        switch (R1C) {
        case ('a') :
        case ('A') :
            for (std::size_t i = 0; i <= qap.numCircuitInputs(); ++i)
                m_vec[3 + i] = (*uit) * T(i + 1);
        }

        for (const auto& constraint : qap.constraintSystem().constraints()) {
            ++uit;

            for (const auto& term : constraint.combo(R1C).terms())
                m_vec[3 + term.index()] += (*uit) * term.coeff();
        }

        for (const auto& v : m_vec)
            if (! v.isZero()) ++m_nonzeroCount;
    }

    std::size_t nonzeroCount() const { return m_nonzeroCount; }
    const std::vector<T>& vec() const { return m_vec; }

    // only used by QAP_QueryIC<T>
    void zeroElement(const std::size_t index) {
        m_vec[index] = T::zero();
    }

private:
    std::size_t m_nonzeroCount;
    std::vector<T> m_vec;
};

template <typename T> using QAP_QueryA = QAP_QueryABC<T, 'A', 0>;
template <typename T> using QAP_QueryB = QAP_QueryABC<T, 'B', 1>;
template <typename T> using QAP_QueryC = QAP_QueryABC<T, 'C', 2>;

////////////////////////////////////////////////////////////////////////////////
// query vector H
//

template <typename T>
class QAP_QueryH
{
public:
    QAP_QueryH(const QAP_SystemPoint<T>& qap)
        : m_nonzeroCount(0),
          m_vec(qap.degree() + 1, T::zero())
    {
        auto ti = T::one();

        for (auto& r : m_vec) {
            r = ti;
            ti *= qap.point();
        }

        for (const auto& v : m_vec)
            if (! v.isZero()) ++m_nonzeroCount;
    }

    std::size_t nonzeroCount() const { return m_nonzeroCount; }
    const std::vector<T>& vec() const { return m_vec; }

private:
    std::size_t m_nonzeroCount;
    std::vector<T> m_vec;
};

////////////////////////////////////////////////////////////////////////////////
// window table dimensions
//

template <typename T>
std::size_t g1_exp_count(const QAP_SystemPoint<T>& qap,
                         const QAP_QueryA<T>& At,
                         const QAP_QueryB<T>& Bt,
                         const QAP_QueryC<T>& Ct,
                         const QAP_QueryH<T>& Ht) {
    return
        2 * (At.nonzeroCount() - qap.numCircuitInputs() + Ct.nonzeroCount())
        + Bt.nonzeroCount()
        + Ht.nonzeroCount()
        + 3 + qap.numVariables() + 1; // K_query.size()
}

template <typename T>
std::size_t g2_exp_count(const QAP_QueryB<T>& Bt) {
    return Bt.nonzeroCount();
}

////////////////////////////////////////////////////////////////////////////////
// randomness derived input consistency coefficients
//

template <typename T>
class QAP_QueryIC
{
public:
    QAP_QueryIC(const QAP_SystemPoint<T>& qap,
                QAP_QueryA<T>& At,
                const T& random_A)
        : m_vec(qap.numCircuitInputs() + 1, T::zero())
    {
        // circuit inputs from At query vector
        for (std::size_t i = 0; i < m_vec.size(); ++i) {
            m_vec[i] = At.vec()[3 + i] * random_A;
#ifdef USE_ASSERT
            assert(! m_vec[i].isZero());
#endif
            At.zeroElement(3 + i);
        }
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
};

////////////////////////////////////////////////////////////////////////////////
// randomness derived vector K
//

template <typename T>
class QAP_QueryK
{
public:
    QAP_QueryK(const QAP_SystemPoint<T>& qap,
               const QAP_QueryA<T>& At,
               const QAP_QueryB<T>& Bt,
               const QAP_QueryC<T>& Ct,
               const T& random_A,
               const T& random_B,
               const T& random_beta)
        : m_vec(3 + qap.numVariables() + 1, T::zero())
    {
        const auto random_C = random_A * random_B;

        for (std::size_t i = 0; i < m_vec.size(); ++i) {
            m_vec[i] = random_beta * (random_A * At.vec()[i]
                                      + random_B * Bt.vec()[i]
                                      + random_C * Ct.vec()[i]);
        }
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
};

////////////////////////////////////////////////////////////////////////////////
// witness vectors A, B, C
//

template <typename T, char R1C>
class QAP_WitnessABC
{
public:
    QAP_WitnessABC(const QAP_SystemPoint<T>& qap,
                   const R1Witness<T>& witness)
        : m_qap(qap),
          m_vec(qap.degree(), T::zero())
    {
        auto uit = m_vec.begin();

        // input consistency
        switch (R1C) {
        case ('a') :
        case ('A') :
            *uit = T::one();
            for (std::size_t i = 0; i < qap.numCircuitInputs(); ++i) {
                *uit += witness[i] * T(i + 2);
            }
        }

        for (const auto& constraint : qap.constraintSystem().constraints()) {
            ++uit;
            *uit += constraint.combo(R1C) * witness;
        }

        qap.FFT()->iFFT(m_vec);
    }

    void cosetFFT() {
        m_qap.FFT()->cosetFFT(m_vec, T::params.multiplicative_generator());
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    const QAP_SystemPoint<T>& m_qap;
    std::vector<T> m_vec;
};

template <typename T> using QAP_WitnessA = QAP_WitnessABC<T, 'A'>;
template <typename T> using QAP_WitnessB = QAP_WitnessABC<T, 'B'>;
template <typename T> using QAP_WitnessC = QAP_WitnessABC<T, 'C'>;

////////////////////////////////////////////////////////////////////////////////
// witness vector H
//

template <typename T>
class QAP_WitnessH
{
public:
    // regular H
    QAP_WitnessH(const QAP_SystemPoint<T>& qap,
                 const QAP_WitnessA<T>& aA, // before cosetFFT()
                 const QAP_WitnessB<T>& aB, // before cosetFFT()
                 const T& random_d1,
                 const T& random_d2,
                 const T& random_d3)
        : m_vec(qap.degree() + 1, T::zero())
    {
        for (std::size_t i = 0; i < qap.degree(); ++i)
            m_vec[i] = random_d2 * aA.vec()[i] + random_d1 * aB.vec()[i];

        m_vec[0] -= random_d3;
        qap.FFT()->add_poly_Z(random_d1 * random_d2, m_vec);
    }

    // temporary H
    QAP_WitnessH(const QAP_SystemPoint<T>& qap,
                 const QAP_WitnessA<T>& aA, // after cosetFFT()
                 const QAP_WitnessB<T>& aB, // after cosetFFT()
                 const QAP_WitnessC<T>& aC) // after cosetFFT()
        : m_vec(qap.degree(), T::zero())
    {
        for (std::size_t i = 0; i < m_vec.size(); ++i)
            m_vec[i] = aA.vec()[i] * aB.vec()[i] - aC.vec()[i];

        qap.FFT()->divide_by_Z_on_coset(m_vec);
        qap.FFT()->icosetFFT(m_vec, T::params.multiplicative_generator());
    }

    void addTemporary(const QAP_WitnessH& tmpH) {
        // make sure to add temporary H, not regular H
#ifdef USE_ASSERT
        assert(tmpH.vec().size() < m_vec.size());
#endif

        for (std::size_t i = 0; i < tmpH.vec().size(); ++i)
            m_vec[i] += tmpH.vec()[i];
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
};

} // namespace snarklib

#endif
