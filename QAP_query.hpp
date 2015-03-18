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

template <template <typename> class SYS, typename T>
class QAP_QueryABC
{
public:
    enum VecSelect { A = 0x001, B = 0x010, C = 0x100 };

    // all vectors A, B, C
    QAP_QueryABC(const QAP_SystemPoint<SYS, T>& qap,
                 const unsigned int mask = A | B | C)
        : m_A(mask & A),
          m_B(mask & B),
          m_C(mask & C),
          m_nonzeroA(0),
          m_nonzeroB(0),
          m_nonzeroC(0),
          m_vecA(m_A ? 3 + qap.numVariables() + 1 : 0, T::zero()),
          m_vecB(m_B ? 3 + qap.numVariables() + 1 : 0, T::zero()),
          m_vecC(m_C ? 3 + qap.numVariables() + 1 : 0, T::zero()),
          m_uit(qap.lagrange_coeffs().begin()),
          m_error(false)
    {
        if (m_A) m_vecA[0] = qap.compute_Z();
        if (m_B) m_vecB[1] = qap.compute_Z();
        if (m_C) m_vecC[2] = qap.compute_Z();

        // input consistency (for A only)
        if (m_A) {
            for (std::size_t i = 0; i <= qap.numCircuitInputs(); ++i)
                m_vecA[3 + i] = (*m_uit) * T(i + 1);
        }

        constraintLoop(qap.constraintSystem());

        if (m_A) {
            m_nonzeroA = std::count_if(m_vecA.begin(),
                                       m_vecA.end(),
                                       [] (const T& v) -> bool {
                                           return ! v.isZero();
                                       });
        }

        if (m_B) {
            m_nonzeroB = std::count_if(m_vecB.begin(),
                                       m_vecB.end(),
                                       [] (const T& v) -> bool {
                                           return ! v.isZero();
                                       });
        }

        if (m_C) {
            m_nonzeroC = std::count_if(m_vecC.begin(),
                                       m_vecC.end(),
                                       [] (const T& v) -> bool {
                                           return ! v.isZero();
                                       });
        }
    }

    std::size_t nonzeroA() const { return m_nonzeroA; }
    std::size_t nonzeroB() const { return m_nonzeroB; }
    std::size_t nonzeroC() const { return m_nonzeroC; }

    std::size_t nonzeroCount() const {
        std::size_t count = 0;

        if (m_A) return count += nonzeroA();
        if (m_B) return count += nonzeroB();
        if (m_C) return count += nonzeroC();

        return count;
    }

    const std::vector<T>& vecA() const { return m_vecA; }
    const std::vector<T>& vecB() const { return m_vecB; }
    const std::vector<T>& vecC() const { return m_vecC; }

    const std::vector<T>& vec() const {
        if (m_A)
            return vecA();
        else if (m_B)
            return vecB();
        else
            return vecC();
    }

    bool operator! () const { return m_error; }

private:
    void constraintLoop(const R1System<T>& S)
    {
        if (m_A && !m_B && !m_C) {
            // only A
            for (const auto& c : S.constraints()) {
                ++m_uit;
                const auto& u = *m_uit;
                for (const auto& t : c.a().terms())
                    m_vecA[3 + t.index()] += u * t.coeff();
            }

        } else if (m_B && !m_A && !m_C) {
            // only B
            for (const auto& c : S.constraints()) {
                ++m_uit;
                const auto& u = *m_uit;
                for (const auto& t : c.b().terms())
                    m_vecB[3 + t.index()] += u * t.coeff();
            }

        } else if (m_C && !m_A && !m_C) {
            // only C
            for (const auto& c : S.constraints()) {
                ++m_uit;
                const auto& u = *m_uit;
                for (const auto& t : c.c().terms())
                    m_vecC[3 + t.index()] += u * t.coeff();
            }

        } else if (m_A && m_B && m_C) {
            // all query vectors
            for (const auto& c : S.constraints()) {
                ++m_uit;
                const auto& u = *m_uit;

                for (const auto& t : c.a().terms())
                    m_vecA[3 + t.index()] += u * t.coeff();

                for (const auto& t : c.b().terms())
                    m_vecB[3 + t.index()] += u * t.coeff();

                for (const auto& t : c.c().terms())
                    m_vecC[3 + t.index()] += u * t.coeff();
            }

        } else {
            // subset of query vectors
            for (const auto& c : S.constraints()) {
                ++m_uit;
                const auto& u = *m_uit;

                if (m_A) {
                    for (const auto& t : c.a().terms())
                        m_vecA[3 + t.index()] += u * t.coeff();
                }

                if (m_B) {
                    for (const auto& t : c.b().terms())
                        m_vecB[3 + t.index()] += u * t.coeff();
                }

                if (m_C) {
                    for (const auto& t : c.c().terms())
                        m_vecC[3 + t.index()] += u * t.coeff();
                }
            }
        }
    }

    void constraintLoop(const HugeSystem<T>& S) {
        m_error = S.mapLambda(
            [this] (const R1System<T>& a) -> bool {
                this->constraintLoop(a);
                return false; // do not write back to disk
            });
    }

    const bool m_A, m_B, m_C;
    std::size_t m_nonzeroA, m_nonzeroB, m_nonzeroC;
    std::vector<T> m_vecA, m_vecB, m_vecC;
    typename std::vector<T>::const_iterator m_uit;
    bool m_error;
};

template <template <typename> class SYS, typename T>
class QAP_QueryA : public QAP_QueryABC<SYS, T>
{
public:
    QAP_QueryA(const QAP_SystemPoint<SYS, T>& qap)
        : QAP_QueryABC<SYS, T>{qap, QAP_QueryABC<SYS, T>::VecSelect::A}
    {}
};

template <template <typename> class SYS, typename T>
class QAP_QueryB : public QAP_QueryABC<SYS, T>
{
public:
    QAP_QueryB(const QAP_SystemPoint<SYS, T>& qap)
        : QAP_QueryABC<SYS, T>{qap, QAP_QueryABC<SYS, T>::VecSelect::B}
    {}
};

template <template <typename> class SYS, typename T>
class QAP_QueryC : public QAP_QueryABC<SYS, T>
{
public:
    QAP_QueryC(const QAP_SystemPoint<SYS, T>& qap)
        : QAP_QueryABC<SYS, T>{qap, QAP_QueryABC<SYS, T>::VecSelect::C}
    {}
};

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
                         const QAP_QueryABC<SYS, T>& ABCt,
                         const QAP_QueryH<SYS, T>& Ht) {
    return g1_exp_count(qap.numVariables(),
                        qap.numCircuitInputs(),
                        ABCt.nonzeroA(),
                        ABCt.nonzeroB(),
                        ABCt.nonzeroC(),
                        Ht.nonzeroCount());
}

template <typename T>
T g2_exp_count(const T Bt_nonzeroCount) {
    return Bt_nonzeroCount;
}

template <template <typename> class SYS, typename T>
std::size_t g2_exp_count(const QAP_QueryABC<SYS, T>& ABCt) {
    return g2_exp_count(ABCt.nonzeroB());
}

////////////////////////////////////////////////////////////////////////////////
// randomness derived input consistency coefficients
//

template <template <typename> class SYS, typename T>
std::vector<T> qap_query_IC(const QAP<T>& qap,
                            const QAP_QueryABC<SYS, T>& ABCt,
                            const T& random_rA)
{
    std::vector<T> vec(qap.numCircuitInputs() + 1, T::zero());

    // circuit inputs from At query vector
    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec[i] = ABCt.vecA()[3 + i] * random_rA;

#ifdef USE_ASSERT
        assert(! vec[i].isZero());
#endif
    }

    return vec;
}

template <template <typename> class SYS, typename T>
std::vector<T> qap_query_IC(const QAP<T>& qap,
                            const QAP_QueryABC<SYS, T>& ABCt)
{
    auto vec = ABCt.vecA();

    // zero out the circuit inputs
    for (std::size_t i = 0; i < qap.numCircuitInputs() + 1; ++i) {
        vec[3 + i] = T::zero();
    }

    return vec;
}

template <typename T>
class QAP_QueryIC
{
public:
    // use with accumVector() to avoid having all vectors in memory
    QAP_QueryIC(const QAP<T>& qap,
                const T& random_rA)
        : m_vec(qap.numCircuitInputs() + 1, T::zero()),
          m_random_rA(random_rA)
    {}

    // blinded random_A is windowed exponentiation table generator
    QAP_QueryIC(const QAP<T>& qap)
        : QAP_QueryIC{qap, T::one()}
    {}

    // only need to accumulate blocks within number of circuit inputs
    // returns true if At is modified, false otherwise
    bool accumVector(BlockVector<T>& At) {
        const auto limit = std::min(At.stopIndex(), m_vec.size());
        if (At.startIndex() >= limit) {
            return false;

        } else {
            for (std::size_t i = At.startIndex(); i < limit; ++i) {
                m_vec[i] = At[3 + i] * m_random_rA;
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
    const T m_random_rA;
};

////////////////////////////////////////////////////////////////////////////////
// randomness derived vector K
//

template <template <typename> class SYS, typename T>
std::vector<T> qap_query_K(const QAP<T>& qap,
                           const QAP_QueryABC<SYS, T>& ABCt,
                           const T& random_beta_rA,
                           const T& random_beta_rB,
                           const T& random_beta_rC)
{
    std::vector<T> vec(3 + qap.numVariables() + 1, T::zero());

    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec[i] =
            random_beta_rA * ABCt.vecA()[i] +
            random_beta_rB * ABCt.vecB()[i] +
            random_beta_rC * ABCt.vecC()[i];
    }

    return vec;
}

template <typename T>
class QAP_QueryK
{
public:
    // use with accumVector() to avoid having all vectors in memory
    QAP_QueryK(const QAP<T>& qap,
               const T& random_beta_rA,
               const T& random_beta_rB,
               const T& random_beta_rC)
        : m_vec(3 + qap.numVariables() + 1, T::zero()),
          m_random_beta_rA(random_beta_rA),
          m_random_beta_rB(random_beta_rB),
          m_random_beta_rC(random_beta_rC)
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
            m_vec[i] =
                m_random_beta_rA * At[i] +
                m_random_beta_rB * Bt[i] +
                m_random_beta_rC * Ct[i];
        }
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
    const T& m_random_beta_rA;
    const T& m_random_beta_rB;
    const T& m_random_beta_rC;
};

} // namespace snarklib

#endif
