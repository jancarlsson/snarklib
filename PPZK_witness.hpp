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
#include "QAP.hpp"
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// witnesses A, B, C
//

template <typename GA, typename GB, typename FR, std::size_t Z_INDEX>
class PPZK_WitnessABC
{
public:
    PPZK_WitnessABC(const QAP_SystemPoint<FR>& qap,
                    const R1Witness<FR>& witness,
                    const FR& random_d)
        : m_numVariables(qap.numVariables()),
          m_witness(*witness),
          m_random_d(random_d)
    {}

    void accumQuery(const SparseVector<Pairing<GA, GB>>& query,
                    const std::size_t reserveTune,
                    ProgressCallback* callback) {
        m_val = m_val
            + m_random_d * query.getElementForIndex(Z_INDEX)
            + query.getElementForIndex(3);

        if (0 == reserveTune) {
            m_val = m_val + multiExp01(query,
                                       m_witness,
                                       4,
                                       4 + m_numVariables,
                                       callback);
        } else {
            m_val = m_val + multiExp01(query,
                                       m_witness,
                                       4,
                                       4 + m_numVariables,
                                       m_numVariables / reserveTune,
                                       callback);
        }
    }

    void accumQuery(const SparseVector<Pairing<GA, GB>>& query,
                    ProgressCallback* callback = nullptr) {
        accumQuery(query, 0, callback);
    }

    const Pairing<GA, GB>& val() const { return m_val; }

private:
    const std::size_t m_numVariables;
    const std::vector<FR>& m_witness;
    const FR& m_random_d;
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
    PPZK_WitnessK(const R1Witness<Fr>& witness,
                  const Fr& random_d1,
                  const Fr& random_d2,
                  const Fr& random_d3)
        : m_witness(*witness),
          m_random_d1(random_d1),
          m_random_d2(random_d2),
          m_random_d3(random_d3)
    {}

    void accumQuery(const BlockVector<G1>& query,
                    const std::size_t reserveTune,
                    ProgressCallback* callback = nullptr)
    {
        if (0 == query.block()[0]) {
#ifdef USE_ASSERT
            assert(query.size() >= 4);
#endif

            m_val = m_val
                + m_random_d1 * query[0]
                + m_random_d2 * query[1]
                + m_random_d3 * query[2]
                + query[3];

            if (0 == reserveTune) {
                m_val = m_val + multiExp01(
                    std::vector<G1>(query.vec().begin() + 4, query.vec().end()),
                    m_witness,
                    callback);

            } else {
                m_val = m_val + multiExp01(
                    std::vector<G1>(query.vec().begin() + 4, query.vec().end()),
                    m_witness,
                    (query.vec().size() - 4) / reserveTune,
                    callback);
            }

        } else {
            if (0 == reserveTune) {
                m_val = m_val + multiExp01(query.vec(),
                                           m_witness,
                                           callback);
            } else {
                m_val = m_val + multiExp01(query.vec(),
                                           m_witness,
                                           query.vec().size() / reserveTune,
                                           callback);
            }
        }
    }

    const G1& val() const { return m_val; }

private:
    const std::vector<Fr>& m_witness;
    const Fr& m_random_d1;
    const Fr& m_random_d2;
    const Fr& m_random_d3;
    G1 m_val;
};

} // namespace snarklib

#endif
