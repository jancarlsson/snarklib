#ifndef _SNARKLIB_PPZK_WITNESS_HPP_
#define _SNARKLIB_PPZK_WITNESS_HPP_

#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
#include "AuxSTL.hpp"
#include "MultiExp.hpp"
#include "Pairing.hpp"
#include "ProgressCallback.hpp"
#include "QAP_system.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// witnesses A, B, C
//

template <typename GA, typename GB, typename FR, std::size_t Z_INDEX>
class PPZK_WitnessABC
{
public:
    typedef SparseVector<Pairing<GA, GB>> SparseVec;
    typedef Pairing<GA, GB> Val;

    PPZK_WitnessABC()
        : m_witness(nullptr),
          m_random_d(nullptr)
    {}

    PPZK_WitnessABC(const std::size_t qapNumVariables,
                    const R1Witness<FR>& witness,
                    const FR& random_d)
        : m_numVariables(qapNumVariables),
          m_witness(std::addressof(*witness)),
          m_random_d(std::addressof(random_d))
    {}

    PPZK_WitnessABC(const QAP<FR>& qap,
                    const R1Witness<FR>& witness,
                    const FR& random_d)
        : PPZK_WitnessABC{qap.numVariables(), witness, random_d}
    {}

    void accumQuery(const SparseVector<Pairing<GA, GB>>& query,
                    const std::size_t reserveTune,
                    ProgressCallback* callback) {
        m_val = m_val
            + (*m_random_d) * query.getElementForIndex(Z_INDEX)
            + query.getElementForIndex(3)
            + multiExp01(query,
                         *m_witness,
                         4,
                         4 + m_numVariables,
                         0 == reserveTune ? reserveTune : m_numVariables / reserveTune,
                         callback);
    }

    void accumQuery(const SparseVector<Pairing<GA, GB>>& query,
                    ProgressCallback* callback = nullptr) {
        accumQuery(query, 0, callback);
    }

    const Pairing<GA, GB>& val() const { return m_val; }

private:
    std::size_t m_numVariables;
    const std::vector<FR>* m_witness;
    const FR* m_random_d;
    Pairing<GA, GB> m_val;
};

template <typename PAIRING> using PPZK_WitnessA =
    PPZK_WitnessABC<typename PAIRING::G1, typename PAIRING::G1, typename PAIRING::Fr, 0>;

template <typename PAIRING> using PPZK_WitnessB =
    PPZK_WitnessABC<typename PAIRING::G2, typename PAIRING::G1, typename PAIRING::Fr, 1>;

template <typename PAIRING> using PPZK_WitnessC =
    PPZK_WitnessABC<typename PAIRING::G1, typename PAIRING::G1, typename PAIRING::Fr, 2>;

////////////////////////////////////////////////////////////////////////////////
// witness H
//

template <typename PAIRING>
class PPZK_WitnessH
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;

public:
    PPZK_WitnessH() = default;

    void accumQuery(const BlockVector<G1>& query,
                    const BlockVector<Fr>& scalar,
                    ProgressCallback* callback = nullptr)
    {
#ifdef USE_ASSERT
        assert(query.space() == scalar.space() &&
               query.block() == scalar.block());
#endif

        m_val = m_val + multiExp(query.vec(),
                                 scalar.vec(),
                                 callback);
    }

    const G1& val() const { return m_val; }

private:
    G1 m_val;
};

////////////////////////////////////////////////////////////////////////////////
// witness K
//

template <typename PAIRING>
class PPZK_WitnessK
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;

public:
    PPZK_WitnessK()
        : m_witness(nullptr),
          m_random_d1(nullptr),
          m_random_d2(nullptr),
          m_random_d3(nullptr)
    {}

    PPZK_WitnessK(const R1Witness<Fr>& witness,
                  const Fr& random_d1,
                  const Fr& random_d2,
                  const Fr& random_d3)
        : m_witness(std::addressof(*witness)),
          m_random_d1(std::addressof(random_d1)),
          m_random_d2(std::addressof(random_d2)),
          m_random_d3(std::addressof(random_d3))
    {}

    void accumQuery(const BlockVector<G1>& query,
                    const std::size_t reserveTune,
                    ProgressCallback* callback = nullptr)
    {
        std::size_t startOffset = 0;

        if (0 == query.block()[0]) {
#ifdef USE_ASSERT
            assert(query.size() >= 4);
#endif

            m_val = m_val
                + (*m_random_d1) * query[0]
                + (*m_random_d2) * query[1]
                + (*m_random_d3) * query[2]
                + query[3];

            startOffset = 4;
        }

        m_val = m_val + multiExp01(
            query,
            startOffset,
            4,
            *m_witness,
            0 == reserveTune ? reserveTune : (query.size() - startOffset) / reserveTune,
            callback);
    }

    const G1& val() const { return m_val; }

private:
    const std::vector<Fr>* m_witness;
    const Fr* m_random_d1;
    const Fr* m_random_d2;
    const Fr* m_random_d3;
    G1 m_val;
};

} // namespace snarklib

#endif
