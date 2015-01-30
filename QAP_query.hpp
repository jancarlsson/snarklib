#ifndef _SNARKLIB_QAP_QUERY_HPP_
#define _SNARKLIB_QAP_QUERY_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
#include "AuxSTL.hpp"
#include "HugeSystem.hpp"
#include "QAP_system.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// query vectors A, B, C
//

template <template <typename> class SYS, typename T, char R1C, std::size_t Z_INDEX>
class QAP_QueryABC
{
public:
    QAP_QueryABC(const QAP_SystemPoint<SYS, T>& qap)
        : m_nonzeroCount(0),
          m_vec(3 + qap.numVariables() + 1, T::zero()),
          m_uit(qap.lagrange_coeffs().begin()),
          m_error(false)
    {
        m_vec[Z_INDEX] = qap.compute_Z();

        // input consistency
        switch (R1C) {
        case ('a') :
        case ('A') : 
            for (std::size_t i = 0; i <= qap.numCircuitInputs(); ++i)
                m_vec[3 + i] = (*m_uit) * T(i + 1);
        }

        constraintLoop(qap.constraintSystem());

        for (const auto& v : m_vec)
            if (! v.isZero()) ++m_nonzeroCount;
    }

    std::size_t nonzeroCount() const { return m_nonzeroCount; }
    const std::vector<T>& vec() const { return m_vec; }

    bool operator! () const { return m_error; }

    // only used by QAP_QueryIC<T>
    void zeroElement(const std::size_t index) {
        m_vec[index] = T::zero();
    }

private:
    void constraintLoop(const R1System<T>& S) {
        for (const auto& constraint : S.constraints()) {
            ++m_uit;

            for (const auto& term : constraint.combo(R1C).terms())
                m_vec[3 + term.index()] += (*m_uit) * term.coeff();
        }
    }

    void constraintLoop(const HugeSystem<T>& S) {
        m_error = S.mapLambda(
            [this] (const R1System<T>& a) -> bool {
                this->constraintLoop(a);
                return false; // do not write back to disk
            });
    }

    std::size_t m_nonzeroCount;
    std::vector<T> m_vec;
    typename std::vector<T>::const_iterator m_uit;
    bool m_error;
};

template <template <typename> class SYS, typename T> using
QAP_QueryA = QAP_QueryABC<SYS, T, 'A', 0>;

template <template <typename> class SYS, typename T> using
QAP_QueryB = QAP_QueryABC<SYS, T, 'B', 1>;

template <template <typename> class SYS, typename T> using
QAP_QueryC = QAP_QueryABC<SYS, T, 'C', 2>;

////////////////////////////////////////////////////////////////////////////////
// query vector H
//

template <template <typename> class SYS, typename T>
class QAP_QueryH
{
public:
    QAP_QueryH(const QAP_SystemPoint<SYS, T>& qap)
        : m_nonzeroCount(0),
          m_vec(qap.degree() + 1, T::zero())
    {
        auto ti = T::one();

        for (auto& r : m_vec) {
            r = ti;

            if (! r.isZero()) ++m_nonzeroCount;

            ti *= qap.point();
        }
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
T g1_exp_count(const T qap_numVariables,
               const T qap_numCircuitInputs,
               const T At_nonzeroCount,
               const T Bt_nonzeroCount,
               const T Ct_nonzeroCount,
               const T Ht_nonzeroCount) {
    return
        2 * (At_nonzeroCount - qap_numCircuitInputs + Ct_nonzeroCount)
        + Bt_nonzeroCount
        + Ht_nonzeroCount
        + 3 + qap_numVariables + 1; // K_query.size()
}

template <template <typename> class SYS, typename T>
std::size_t g1_exp_count(const QAP_SystemPoint<SYS, T>& qap,
                         const QAP_QueryA<SYS, T>& At,
                         const QAP_QueryB<SYS, T>& Bt,
                         const QAP_QueryC<SYS, T>& Ct,
                         const QAP_QueryH<SYS, T>& Ht) {
    return g1_exp_count(qap.numVariables(),
                        qap.numCircuitInputs(),
                        At.nonzeroCount(),
                        Bt.nonzeroCount(),
                        Ct.nonzeroCount(),
                        Ht.nonzeroCount());
}

template <typename T>
T g2_exp_count(const T Bt_nonzeroCount) {
    return Bt_nonzeroCount;
}

template <template <typename> class SYS, typename T>
std::size_t g2_exp_count(const QAP_QueryB<SYS, T>& Bt) {
    return g2_exp_count(Bt.nonzeroCount());
}

////////////////////////////////////////////////////////////////////////////////
// randomness derived input consistency coefficients
//

template <template <typename> class SYS, typename T>
class QAP_QueryIC
{
public:
    QAP_QueryIC(const QAP<T>& qap,
                QAP_QueryA<SYS, T>& At,
                const T& random_A)
        : m_vec(qap.numCircuitInputs() + 1, T::zero()),
          m_random_A(random_A)
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

    // use with accumVector() to avoid having all vectors in memory
    QAP_QueryIC(const QAP<T>& qap,
                const T& random_A)
        : m_vec(qap.numCircuitInputs() + 1, T::zero()),
          m_random_A(random_A)
    {}

    // only need to accumulate blocks within number of circuit inputs
    // returns true if At is modified, false otherwise
    bool accumVector(BlockVector<T>& At) {
        const auto limit = std::min(At.stopIndex(), m_vec.size());
        if (At.startIndex() >= limit) {
            return false;

        } else {
            for (std::size_t i = At.startIndex(); i < limit; ++i) {
                m_vec[i] = At[3 + i] * m_random_A;
#ifdef USE_ASSERT
                assert(! m_vec[i].isZero());
#endif
                At[3 + i] = T::zero();
            }

            return true;
        }
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
    const T& m_random_A;
};

////////////////////////////////////////////////////////////////////////////////
// randomness derived vector K
//

template <template <typename> class SYS, typename T>
class QAP_QueryK
{
public:
    QAP_QueryK(const QAP<T>& qap,
               const QAP_QueryA<SYS, T>& At,
               const QAP_QueryB<SYS, T>& Bt,
               const QAP_QueryC<SYS, T>& Ct,
               const T& random_A,
               const T& random_B,
               const T& random_beta)
        : m_vec(3 + qap.numVariables() + 1, T::zero()),
          m_random_A(random_A),
          m_random_B(random_B),
          m_random_C(random_A * random_B),
          m_random_beta(random_beta)
    {
        for (std::size_t i = 0; i < m_vec.size(); ++i) {
            m_vec[i] = random_beta * (random_A * At.vec()[i]
                                      + random_B * Bt.vec()[i]
                                      + m_random_C * Ct.vec()[i]);
        }
    }

    // use with accumVector() to avoid having all vectors in memory
    QAP_QueryK(const QAP<T>& qap,
               const T& random_A,
               const T& random_B,
               const T& random_beta)
        : m_vec(3 + qap.numVariables() + 1, T::zero()),
          m_random_A(random_A),
          m_random_B(random_B),
          m_random_C(random_A * random_B),
          m_random_beta(random_beta)
    {}

    // must accumulate all blocks
    void accumVector(const BlockVector<T>& At,
                     const BlockVector<T>& Bt,
                     const BlockVector<T>& Ct) {
#ifdef USE_ASSERT
        assert(At.space() == Bt.space() &&
               Bt.space() == Ct.space() &&
               At.block() == Bt.block() &&
               Bt.block() == Ct.block());
#endif

        for (std::size_t i = At.startIndex(); i < At.stopIndex(); ++i) {
            m_vec[i] = m_random_beta * (m_random_A * At[i]
                                        + m_random_B * Bt[i]
                                        + m_random_C * Ct[i]);
        }
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
    const T& m_random_A;
    const T& m_random_B;
    const T m_random_C;
    const T& m_random_beta;
};

} // namespace snarklib

#endif
