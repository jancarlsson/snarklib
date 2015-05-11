#ifndef _SNARKLIB_QAP_SYSTEM_HPP_
#define _SNARKLIB_QAP_SYSTEM_HPP_

#include <cstdint>
#include <vector>

#include <snarklib/HugeSystem.hpp>
#include <snarklib/LagrangeFFT.hpp>
#include <snarklib/LagrangeFFTX.hpp>
#include <snarklib/Rank1DSL.hpp>

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

    QAP(const HugeSystem<T>& hugeSystem,
        const std::size_t numCircuitInputs)
        : m_degree(LagrangeFFT<T>::getDegree(hugeSystem.totalConstraints() + 1)),
          m_numVariables(hugeSystem.maxIndex()),
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

// SYS may be R1System<T> or HugeSystem<T>
template <template <typename> class SYS, typename T>
class QAP_SystemPoint : public QAP<T>
{
public:
    // keypair generation: ABCH, IC coefficients, and K
    QAP_SystemPoint(const SYS<T>& constraintSystem,
                    const std::size_t numCircuitInputs,
                    const T& point)
        : QAP<T>(constraintSystem, numCircuitInputs),
          m_constraintSystem(constraintSystem),
          m_point(point),
          m_compute_Z(QAP<T>::FFT()->compute_Z(point)),
          m_lagrange_coeffs(QAP<T>::FFT()->lagrange_coeffs(point))
    {}

    // proof generation
    QAP_SystemPoint(const SYS<T>& constraintSystem,
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

    const SYS<T>& constraintSystem() const { return m_constraintSystem; }

    const T& point() const { return m_point; }
    const T& compute_Z() const { return m_compute_Z; }
    const std::vector<T>& lagrange_coeffs() const { return m_lagrange_coeffs; }

private:
    const SYS<T>& m_constraintSystem;
    const T m_point, m_compute_Z;
    const std::vector<T> m_lagrange_coeffs;
};

} // namespace snarklib

#endif
