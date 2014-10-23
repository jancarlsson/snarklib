#ifndef _SNARKLIB_QAP_HPP_
#define _SNARKLIB_QAP_HPP_

#include <cassert>
#include <cstdint>
#include <vector>
#include "LagrangeFFT.hpp"
#include "LagrangeFFTX.hpp"
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
// ABCH evaluated at t (from constraint system to generated keypair)
//

template <typename T>
class QAP_ABCH : public QAP<T>
{
public:
    std::size_t degree() const { return QAP<T>::degree(); }
    std::size_t numVariables() const { return QAP<T>::numVariables(); }
    std::size_t numCircuitInputs() const { return QAP<T>::numCircuitInputs(); }
    const LagrangeFFT<T>& FFT() const { return QAP<T>::FFT(); }

    QAP_ABCH(const R1System<T>& constraintSystem,
             const std::size_t numInputs,
             const T& point)
        : QAP<T>(constraintSystem, numInputs),
          m_nonzeroAt(0),
          m_nonzeroBt(0),
          m_nonzeroCt(0),
          m_nonzeroHt(0),
          m_At(3 + numVariables() + 1, T::zero()),
          m_Bt(3 + numVariables() + 1, T::zero()),
          m_Ct(3 + numVariables() + 1, T::zero())
    {
        m_Ht.reserve(degree() + 1);

        const T Z = FFT()->compute_Z(point);

        m_At[0] = Z;
        m_Bt[1] = Z;
        m_Ct[2] = Z;

        const auto u = FFT()->lagrange_coeffs(point);
        auto uit = u.begin();

        for (std::size_t i = 0; i <= numCircuitInputs(); ++i) {
            m_At[3 + i] += (*uit) * T(i + 1);
        }

        for (const auto& constraint : constraintSystem.constraints()) {
            ++uit;

            for (const auto& term : constraint.a().terms()) {
                m_At[3 + term.index()] += (*uit) * term.coeff();
            }

            for (const auto& term : constraint.b().terms()) {
                m_Bt[3 + term.index()] += (*uit) * term.coeff();
            }
            for (const auto& term : constraint.c().terms()) {
                m_Ct[3 + term.index()] += (*uit) * term.coeff();
            }
        }

        auto ti = T::one();
        for (std::size_t i = 0; i <= degree(); ++i) {
            m_Ht.emplace_back(ti);
            ti *= point;
        }

        for (const auto& v : m_At) {
            if (! v.isZero())
                ++m_nonzeroAt;
        }

        for (const auto& v : m_Bt) {
            if (! v.isZero())
                ++m_nonzeroBt;
        }

        for (const auto& v : m_Ct) {
            if (! v.isZero())
                ++m_nonzeroCt;
        }

        for (const auto& v : m_Ht) {
            if (! v.isZero())
                ++m_nonzeroHt;
        }
    }

    const std::vector<T>& A_query() const { return m_At; }
    const std::vector<T>& B_query() const { return m_Bt; }
    const std::vector<T>& C_query() const { return m_Ct; }
    const std::vector<T>& H_query() const { return m_Ht; }

    std::size_t nonzeroAt() const { return m_nonzeroAt; }
    std::size_t nonzeroBt() const { return m_nonzeroBt; }
    std::size_t nonzeroCt() const { return m_nonzeroCt; }
    std::size_t nonzeroHt() const { return m_nonzeroHt; }

    std::vector<T> K_query(const T& beta,
                           const T& rA,
                           const T& rB,
                           const T& rC) const {
        std::vector<T> vec;
        vec.reserve(3 + numVariables() + 1);

        for (std::size_t i = 0; i < 3 + numVariables() + 1; ++i) {
            vec.emplace_back(
                beta * (
                    (rA * A_query()[i]) +
                    (rB * B_query()[i]) +
                    (rC * C_query()[i])));
        }

        return vec;
    }

    std::vector<T> IC_coefficients(const T& rA) {
        std::vector<T> vec;
        vec.reserve(numCircuitInputs() + 1);

        for (std::size_t i = 0; i < numCircuitInputs() + 1; ++i) {
            vec.emplace_back(m_At[3 + i] * rA);
            assert(! vec[i].isZero());
            m_At[3 + i] = T::zero();
        }

        return vec;
    }

    std::size_t g1_exp_count() const {
        return 2 * (nonzeroAt() -
                    numCircuitInputs() +
                    nonzeroCt())
            + nonzeroBt()
            + nonzeroHt()
            + 3 + numVariables() + 1; // K_query.size()
    }
    
    std::size_t g2_exp_count() const {
        return nonzeroBt();
    }

private:
    std::size_t m_nonzeroAt, m_nonzeroBt, m_nonzeroCt, m_nonzeroHt;
    std::vector<T> m_At, m_Bt, m_Ct, m_Ht;
};

////////////////////////////////////////////////////////////////////////////////
// witness map (combined with proving key to construct the proof)
//

template <typename T>
class QAP_Witness : public QAP<T>
{
public:
    std::size_t degree() const { return QAP<T>::degree(); }
    std::size_t numVariables() const { return QAP<T>::numVariables(); }
    std::size_t numCircuitInputs() const { return QAP<T>::numCircuitInputs(); }
    const LagrangeFFT<T>& FFT() const { return QAP<T>::FFT(); }

    QAP_Witness(const R1System<T>& constraintSystem,
                const std::size_t numInputs,
                const std::vector<T>& witness,
                const T& d1,
                const T& d2,
                const T& d3)
        : QAP<T>(constraintSystem, numInputs),
          m_H(degree() + 1, T::zero())
    {
        std::vector<T>
            aA(degree(), T::zero()),
            aB(degree(), T::zero());

        auto
            aAit = aA.begin(),
            aBit = aB.begin();

        *aAit = T::one();
        for (std::size_t i = 0; i < numCircuitInputs(); ++i) {
            *aAit += witness[i] * T(i + 2);
        }

        for (const auto& constraint : constraintSystem.constraints()) {
            ++aAit;
            *aAit += constraint.a() * witness;

            ++aBit;
            *aBit += constraint.b() * witness;
        }

        FFT()->iFFT(aA);
        FFT()->iFFT(aB);

        for (std::size_t i = 0; i < degree(); ++i) {
            m_H[i] = d2 * aA[i] + d1 * aB[i];
        }
        m_H[0] -= d3;
        FFT()->add_poly_Z(d1 * d2, m_H);

        FFT()->cosetFFT(aA, T::params.multiplicative_generator());
        FFT()->cosetFFT(aB, T::params.multiplicative_generator());

        auto& H_tmp = aA;

        for (std::size_t i = 0; i < degree(); ++i) {
            H_tmp[i] = aA[i] * aB[i];
        }

        std::vector<T> aC(degree(), T::zero());
        auto aCit = aC.begin();
        for (const auto& constraint : constraintSystem.constraints()) {
            ++aCit;
            *aCit += constraint.c() * witness;
        }

        FFT()->iFFT(aC);
        FFT()->cosetFFT(aC, T::params.multiplicative_generator());

        for (std::size_t i = 0; i < degree(); ++i) {
            H_tmp[i] = H_tmp[i] - aC[i];
        }

        FFT()->divide_by_Z_on_coset(H_tmp);
        FFT()->icosetFFT(H_tmp, T::params.multiplicative_generator());

        for (std::size_t i = 0; i < degree(); ++i) {
            m_H[i] += H_tmp[i];
        }
    }

    const std::vector<T>& H() const {
        return m_H;
    }

private:
    std::vector<T> m_H;
};

} // namespace snarklib

#endif
