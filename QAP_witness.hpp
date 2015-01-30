#ifndef _SNARKLIB_QAP_WITNESS_HPP_
#define _SNARKLIB_QAP_WITNESS_HPP_

#include <cassert>
#include <cstdint>
#include <vector>
#include "HugeSystem.hpp"
#include "QAP_system.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// witness vectors A, B, C
//

template <template <typename> class SYS, typename T, char R1C>
class QAP_WitnessABC
{
public:
    QAP_WitnessABC(const QAP_SystemPoint<SYS, T>& qap,
                   const R1Witness<T>& witness)
        : m_qap(qap),
          m_witness(witness),
          m_vec(qap.degree(), T::zero()),
          m_uit(m_vec.begin()),
          m_error(false)
    {
        // input consistency
        switch (R1C) {
        case ('a') :
        case ('A') :
            *m_uit = T::one();
            for (std::size_t i = 0; i < qap.numCircuitInputs(); ++i) {
                *m_uit += witness[i] * T(i + 2);
            }
        }

        constraintLoop(qap.constraintSystem());

        qap.FFT()->iFFT(m_vec);
    }

    void cosetFFT() {
        m_qap.FFT()->cosetFFT(m_vec, T::params.multiplicative_generator());
    }

    const std::vector<T>& vec() const { return m_vec; }

    bool operator! () const { return m_error; }

private:
    void constraintLoop(const R1System<T>& S) {
        for (const auto& constraint : S.constraints()) {
            ++m_uit;
            *m_uit += constraint.combo(R1C) * m_witness;
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
    std::vector<T> m_vec;
    typename std::vector<T>::iterator m_uit;
    bool m_error;
};

template <template <typename> class SYS, typename T> using
QAP_WitnessA = QAP_WitnessABC<SYS, T, 'A'>;

template <template <typename> class SYS, typename T> using
QAP_WitnessB = QAP_WitnessABC<SYS, T, 'B'>;

template <template <typename> class SYS, typename T> using
QAP_WitnessC = QAP_WitnessABC<SYS, T, 'C'>;

////////////////////////////////////////////////////////////////////////////////
// witness vector H
//

template <template <typename> class SYS, typename T>
class QAP_WitnessH
{
public:
    ////////////////////////////////////////
    // regular H
    //

    QAP_WitnessH(const QAP_SystemPoint<SYS, T>& qap,
                 const QAP_WitnessA<SYS, T>& aA, // before cosetFFT()
                 const QAP_WitnessB<SYS, T>& aB, // before cosetFFT()
                 const T& random_d1,
                 const T& random_d2,
                 const T& random_d3)
        : QAP_WitnessH{qap, random_d1, random_d2, random_d3}
    {
        for (std::size_t i = 0; i < qap.degree(); ++i)
            m_vec[i] = random_d2 * aA.vec()[i] + random_d1 * aB.vec()[i];

        m_vec[0] -= random_d3;
        qap.FFT()->add_poly_Z(random_d1 * random_d2, m_vec);
    }

    // use with accumVector() to avoid having all vectors in memory
    QAP_WitnessH(const QAP_SystemPoint<SYS, T>& qap,
                 const T& random_d1,
                 const T& random_d2,
                 const T& random_d3)
        : m_vec(qap.degree() + 1, T::zero()),
          m_blockCount(0),
          m_qap(qap),
          m_random_d1(random_d1),
          m_random_d2(random_d2),
          m_random_d3(random_d3)
    {}

    // must accumulate all blocks
    void accumVector(const BlockVector<T>& aA,   // before cosetFFT()
                     const BlockVector<T>& aB) { // before cosetFFT()
#ifdef USE_ASSERT
        assert(aA.space() == aB.space() &&
               aA.block() == aB.block());
#endif

        for (std::size_t i = aA.startIndex(); i < aA.stopIndex(); ++i) {
            m_vec[i] = m_random_d2 * aA[i] + m_random_d1 * aB[i];
        }

        if (aA.space().blockID()[0] == ++m_blockCount) {
            m_vec[0] -= m_random_d3;
            m_qap.FFT()->add_poly_Z(m_random_d1 * m_random_d2, m_vec);
        }
    }

    ////////////////////////////////////////
    // temporary H
    //

    QAP_WitnessH(const QAP_SystemPoint<SYS, T>& qap,
                 const QAP_WitnessA<SYS, T>& aA, // after cosetFFT()
                 const QAP_WitnessB<SYS, T>& aB, // after cosetFFT()
                 const QAP_WitnessC<SYS, T>& aC) // after cosetFFT()
        : QAP_WitnessH{qap}
    {
        for (std::size_t i = 0; i < m_vec.size(); ++i)
            m_vec[i] = aA.vec()[i] * aB.vec()[i] - aC.vec()[i];

        qap.FFT()->divide_by_Z_on_coset(m_vec);
        qap.FFT()->icosetFFT(m_vec, T::params.multiplicative_generator());
    }

    // use with accumVector() to avoid having all vectors in memory
    QAP_WitnessH(const QAP_SystemPoint<SYS, T>& qap)
        : m_vec(qap.degree(), T::zero()),
          m_blockCount(0),
          m_qap(qap),
          m_random_d1(T::zero()), // not used
          m_random_d2(T::zero()), // not used
          m_random_d3(T::zero())  // not used
    {}

    // must accumulate all blocks
    void accumVector(const BlockVector<T>& aA,   // after cosetFFT()
                     const BlockVector<T>& aB,   // after cosetFFT()
                     const BlockVector<T>& aC) { // after cosetFFT()
#ifdef USE_ASSERT
        assert(aA.space() == aB.space() &&
               aB.space() == aC.space() &&
               aA.block() == aB.block() &&
               aB.block() == aC.block());
#endif

        for (std::size_t i = aA.startIndex(); i < aA.stopIndex(); ++i) {
            m_vec[i] = aA[i] * aB[i] - aC[i];
        }

        if (aA.space().blockID()[0] == ++m_blockCount) {
            m_qap.FFT()->divide_by_Z_on_coset(m_vec);
            m_qap.FFT()->icosetFFT(m_vec, T::params.multiplicative_generator());
        }
    }

    // add regular and temporary H together
    void addTemporary(const QAP_WitnessH& tmpH) {
#ifdef USE_ASSERT
        // make sure to add temporary H, not regular H
        assert(tmpH.vec().size() < m_vec.size());
#endif

        for (std::size_t i = 0; i < tmpH.vec().size(); ++i)
            m_vec[i] += tmpH.vec()[i];
    }

    const std::vector<T>& vec() const { return m_vec; }

private:
    std::vector<T> m_vec;
    std::size_t m_blockCount;
    const QAP_SystemPoint<SYS, T>& m_qap;
    const T m_random_d1, m_random_d2, m_random_d3;
};

} // namespace snarklib

#endif
