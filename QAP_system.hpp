#ifndef _SNARKLIB_QAP_SYSTEM_HPP_
#define _SNARKLIB_QAP_SYSTEM_HPP_

#include <cstdint>
#include <vector>

#include <snarklib/LagrangeFFTX.hpp>

#ifndef DISABLE_PARNO_SOUNDNESS_FIX
#define PARNO_SOUNDNESS_FIX
#endif

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// base class
//

template <template <typename> class SYS, typename T>
class QAP
{
public:
    virtual ~QAP() = default;

    std::size_t numVariables() const { return m_numVariables; }
    std::size_t numConstraints() const { return m_numConstraints; }
    std::size_t numCircuitInputs() const { return m_numCircuitInputs; }

    std::size_t degree() const { return m_degree; }

    const LagrangeFFT<T>& FFT() const { return m_FFT; }

protected:
    QAP(const SYS<T>& constraintSystem,
        const std::size_t numCircuitInputs)
        : m_numVariables(constraintSystem.maxIndex()),
          m_numConstraints(constraintSystem.size()),
          m_numCircuitInputs(numCircuitInputs),
          m_degree(LagrangeFFT<T>::getDegree(m_numConstraints
#ifdef PARNO_SOUNDNESS_FIX
                                             + numCircuitInputs
#endif
                                             + 1)),
          m_FFT(m_degree)
    {}

private:
    std::size_t m_numVariables;
    std::size_t m_numConstraints;
    std::size_t m_numCircuitInputs;
    std::size_t m_degree;
    LagrangeFFT<T> m_FFT;
};

////////////////////////////////////////////////////////////////////////////////
// constraint system evaluated at point
//

// SYS may be R1System<T> or HugeSystem<T>
template <template <typename> class SYS, typename T>
class QAP_SystemPoint : public QAP<SYS, T>
{
public:
    // keypair generation: ABCH, IC coefficients, and K
    QAP_SystemPoint(const SYS<T>& constraintSystem,
                    const std::size_t numCircuitInputs,
                    const T& point)
        : QAP<SYS, T>(constraintSystem, numCircuitInputs),
          m_weakPoint(false),
          m_constraintSystem(constraintSystem),
          m_point(point),
          m_compute_Z(QAP<SYS, T>::FFT()->compute_Z(point)),
          m_lagrange_coeffs(QAP<SYS, T>::FFT()->lagrange_coeffs(point, m_weakPoint))
    {}

    // proof generation
    QAP_SystemPoint(const SYS<T>& constraintSystem,
                    const std::size_t numCircuitInputs)
        : QAP<SYS, T>(constraintSystem, numCircuitInputs),
          m_weakPoint(false),
          m_constraintSystem(constraintSystem),
          m_point(T::zero()),
          m_compute_Z(T::zero()),
          m_lagrange_coeffs()
    {}

    std::size_t numVariables() const { return QAP<SYS, T>::numVariables(); }
    std::size_t numConstraints() const { return QAP<SYS, T>::numConstraints(); }
    std::size_t numCircuitInputs() const { return QAP<SYS, T>::numCircuitInputs(); }

    std::size_t degree() const { return QAP<SYS, T>::degree(); }

    const SYS<T>& constraintSystem() const { return m_constraintSystem; }

    bool weakPoint() const { return m_weakPoint; }
    const T& point() const { return m_point; }
    const T& compute_Z() const { return m_compute_Z; }
    const std::vector<T>& lagrange_coeffs() const { return m_lagrange_coeffs; }

private:
    bool m_weakPoint;
    const SYS<T>& m_constraintSystem;
    const T m_point, m_compute_Z;
    const std::vector<T> m_lagrange_coeffs;
};

} // namespace snarklib

#endif
