#ifndef _SNARKLIB_PPZK_QUERY_HPP_
#define _SNARKLIB_PPZK_QUERY_HPP_

#include <cstdint>
#include <istream>
#include <ostream>
#include <vector>
#include "AuxSTL.hpp"
#include "Group.hpp"
#include "MultiExp.hpp"
#include "Pairing.hpp"
#include "ProgressCallback.hpp"
#include "Rank1DSL.hpp"
#include "WindowExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// input consistency query vector
//

template <typename PAIRING>
class PPZK_QueryIC
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;

public:
    PPZK_QueryIC() = default;

    // only used for roundtrip marshalling tests
    PPZK_QueryIC(const G1& base,
                 const std::vector<G1>& encoded_terms)
        : m_base(base),
          m_coeffs(),
          m_encoded_terms(encoded_terms)
    {}

    PPZK_QueryIC(const std::vector<Fr>& qap_query) // QAP query IC
        : m_base(qap_query[0] * G1::one()),
          m_coeffs(qap_query.begin() + 1, qap_query.end()),
          m_encoded_terms(qap_query.size() - 1, G1::zero())
    {}

    void accumTable(const WindowExp<G1>& g1_table,
                    ProgressCallback* callback = nullptr) {
        g1_table.batchExp(m_encoded_terms,
                          m_coeffs,
                          callback);
    }

    PPZK_QueryIC accumWitness(const R1Witness<Fr>& witness) const {
        G1 base = m_base;
        std::vector<G1> encoded_terms;

        const std::size_t
            wsize = witness.size(),
            tsize = input_size();

        if (wsize < tsize) {
            base = base + multiExp(
                std::vector<G1>(m_encoded_terms.begin(),
                                m_encoded_terms.begin() + wsize),
                *witness);

            encoded_terms = std::vector<G1>(m_encoded_terms.begin() + wsize,
                                            m_encoded_terms.end());

        } else if (wsize > tsize) {
            base = base + multiExp(m_encoded_terms,
                                   *witness.truncate(tsize));

        } else {
            base = base + multiExp(m_encoded_terms,
                                   *witness);
        }

        return PPZK_QueryIC(base, encoded_terms);
    }

    const G1& base() const {
        return m_base;
    }

    std::size_t input_size() const {
        return m_encoded_terms.size();
    }

    const std::vector<G1>& encoded_terms() const {
        return m_encoded_terms;
    }

    // only used for roundtrip marshalling tests
    bool operator== (const PPZK_QueryIC& other) const {
        return
            base() == other.base() &&
            encoded_terms() == other.encoded_terms();
    }

    // only used for roundtrip marshalling tests
    bool operator!= (const PPZK_QueryIC& other) const {
        return ! (*this == other);
    }

    void marshal_out(std::ostream& os) const {
        base().marshal_out(os);
        snarklib::marshal_out(os, encoded_terms());
    }

    bool marshal_in(std::istream& is) {
        return
            m_base.marshal_in(is) &&
            snarklib::marshal_in(is, m_encoded_terms);
    }

    void clear() {
        m_base = G1::zero();
        m_encoded_terms.clear();
    }

    bool empty() const {
        return
            base().isZero() ||
            encoded_terms().empty();
    }

private:
    G1 m_base;
    std::vector<Fr> m_coeffs;
    std::vector<G1> m_encoded_terms;
};

////////////////////////////////////////////////////////////////////////////////
// query vectors A, B, C
//

template <typename GA, typename GB, typename FR>
class PPZK_QueryABC
{
public:
    PPZK_QueryABC(const BlockVector<FR>& qap_query, // QAP query A, B, C
                  const FR& random_v,
                  const FR& random_alpha)
        : m_qap_query(qap_query),
          m_random_v(random_v),
          m_random_prod(random_v * random_alpha)
    {}

    void accumTable(const WindowExp<GA>& ga_table,
                    const WindowExp<GB>& gb_table,
                    ProgressCallback* callback = nullptr) {
        if (m_vec.empty()) {
            m_vec = batchExp(ga_table,
                             gb_table,
                             m_random_v,
                             m_random_prod,
                             m_qap_query,
                             callback);
        } else {
            batchExp(m_vec,
                     ga_table,
                     gb_table,
                     m_random_v,
                     m_random_prod,
                     m_qap_query,
                     callback);
        }
    }

    const SparseVector<Pairing<GA, GB>>& vec() const { return m_vec; }

private:
    const BlockVector<FR>& m_qap_query;
    const FR& m_random_v;
    const FR m_random_prod;
    SparseVector<Pairing<GA, GB>> m_vec;
};

template <typename PAIRING> using PPZK_QueryA =
    PPZK_QueryABC<typename PAIRING::G1, typename PAIRING::G1, typename PAIRING::Fr>;

template <typename PAIRING> using PPZK_QueryB =
    PPZK_QueryABC<typename PAIRING::G2, typename PAIRING::G1, typename PAIRING::Fr>;

template <typename PAIRING> using PPZK_QueryC =
    PPZK_QueryABC<typename PAIRING::G1, typename PAIRING::G1, typename PAIRING::Fr>;

////////////////////////////////////////////////////////////////////////////////
// query vectors H, K
//

template <typename PAIRING>
class PPZK_QueryHK
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;

public:
    PPZK_QueryHK(const BlockVector<Fr>& qap_query) // QAP query H, K
        : m_qap_query(qap_query),
          m_vec(qap_query.space(), qap_query.block())
    {}

    void accumTable(const WindowExp<G1>& g1_table,
                    ProgressCallback* callback = nullptr) {
        g1_table.batchExp(m_vec,
                          m_qap_query,
                          callback);
    }

    const BlockVector<G1>& vec() const { return m_vec; }
    const std::vector<G1>& vvec() const { return m_vec.vec(); }
    BlockVector<G1>& lvec() { return m_vec; }
    std::vector<G1>& llvec() { return m_vec.lvec(); }

private:
    const BlockVector<Fr>& m_qap_query;
    BlockVector<G1> m_vec;
};

template <typename PAIRING> using PPZK_QueryH = PPZK_QueryHK<PAIRING>;
template <typename PAIRING> using PPZK_QueryK = PPZK_QueryHK<PAIRING>;

} // namespace snarklib

#endif
