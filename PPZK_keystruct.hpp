#ifndef _SNARKLIB_PPZK_KEYSTRUCT_HPP_
#define _SNARKLIB_PPZK_KEYSTRUCT_HPP_

#include <functional>
#include <istream>
#include <ostream>
#include <utility>
#include <vector>

#include <snarklib/AuxSTL.hpp>
#include <snarklib/Group.hpp>
#include <snarklib/Pairing.hpp>
#include <snarklib/PPZK_query.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Proving key
//

template <typename PAIRING>
class PPZK_ProvingKey
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    PPZK_ProvingKey() = default;

    // copy semantics
    PPZK_ProvingKey(const SparseVector<Pairing<G1, G1>>& A_query,
                    const SparseVector<Pairing<G2, G1>>& B_query,
                    const SparseVector<Pairing<G1, G1>>& C_query,
                    const std::vector<G1>& H_query,
                    const std::vector<G1>& K_query)
        : m_A_query(A_query),
          m_B_query(B_query),
          m_C_query(C_query),
          m_H_query(H_query),
          m_K_query(K_query)
    {
        batchSpecial(m_A_query);
        batchSpecial(m_B_query);
        batchSpecial(m_C_query);
        batchSpecial(m_H_query);
        batchSpecial(m_K_query);
    }

    // move semantics
    PPZK_ProvingKey(SparseVector<Pairing<G1, G1>>&& A_query,
                    SparseVector<Pairing<G2, G1>>&& B_query,
                    SparseVector<Pairing<G1, G1>>&& C_query,
                    std::vector<G1>&& H_query,
                    std::vector<G1>&& K_query)
        : m_A_query(std::move(A_query)),
          m_B_query(std::move(B_query)),
          m_C_query(std::move(C_query)),
          m_H_query(std::move(H_query)),
          m_K_query(std::move(K_query))
    {
        batchSpecial(m_A_query);
        batchSpecial(m_B_query);
        batchSpecial(m_C_query);
        batchSpecial(m_H_query);
        batchSpecial(m_K_query);
    }

    const SparseVector<Pairing<G1, G1>>& A_query() const { return m_A_query; }
    const SparseVector<Pairing<G2, G1>>& B_query() const { return m_B_query; }
    const SparseVector<Pairing<G1, G1>>& C_query() const { return m_C_query; }
    const std::vector<G1>& H_query() const { return m_H_query; }
    const std::vector<G1>& K_query() const { return m_K_query; }

    bool operator== (const PPZK_ProvingKey& other) const {
        return
            A_query() == other.A_query() &&
            B_query() == other.B_query() &&
            C_query() == other.C_query() &&
            H_query() == other.H_query() &&
            K_query() == other.K_query();
    }

    bool operator!= (const PPZK_ProvingKey& other) const {
        return ! (*this == other);
    }

    void marshal_out(std::ostream& os) const {
        A_query().marshal_out(os);
        B_query().marshal_out(os);
        C_query().marshal_out(os);
        snarklib::marshal_out(os, H_query());
        snarklib::marshal_out(os, K_query());
    }

    bool marshal_in(std::istream& is) {
        return
            m_A_query.marshal_in(is) &&
            m_B_query.marshal_in(is) &&
            m_C_query.marshal_in(is) &&
            snarklib::marshal_in(is, m_H_query) &&
            snarklib::marshal_in(is, m_K_query);
    }

    void marshal_out_raw(std::ostream& os) const {
        A_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_raw(o);
            });
        
        B_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G2, G1>& a) {
                a.marshal_out_raw(o);
            });

        C_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_raw(o);
            });

        snarklib::marshal_out_raw(os, m_H_query);

        snarklib::marshal_out_raw(os, m_K_query);
    }

    bool marshal_in_raw(std::istream& is) {
        return
            m_A_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_raw(i);
                })
            &&
            m_B_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G2, G1>& a) {
                    return a.marshal_in_raw(i);
                })
            &&
            m_C_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_raw(i);
                })
            &&
            snarklib::marshal_in_raw(is, m_H_query)
            &&
            snarklib::marshal_in_raw(is, m_K_query);
    }

    void marshal_out_special(std::ostream& os) const {
        A_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_special(o);
            });
        
        B_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G2, G1>& a) {
                a.marshal_out_special(o);
            });

        C_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_special(o);
            });

        snarklib::marshal_out_special(os, m_H_query);

        snarklib::marshal_out_special(os, m_K_query);
    }

    bool marshal_in_special(std::istream& is) {
        return
            m_A_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_special(i);
                })
            &&
            m_B_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G2, G1>& a) {
                    return a.marshal_in_special(i);
                })
            &&
            m_C_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_special(i);
                })
            &&
            snarklib::marshal_in_special(is, m_H_query)
            &&
            snarklib::marshal_in_special(is, m_K_query);
    }

    void marshal_out_rawspecial(std::ostream& os) const {
        A_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_rawspecial(o);
            });
        
        B_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G2, G1>& a) {
                a.marshal_out_rawspecial(o);
            });

        C_query().marshal_out(
            os,
            [] (std::ostream& o, const Pairing<G1, G1>& a) {
                a.marshal_out_rawspecial(o);
            });

        snarklib::marshal_out_rawspecial(os, m_H_query);

        snarklib::marshal_out_rawspecial(os, m_K_query);
    }

    bool marshal_in_rawspecial(std::istream& is) {
        return
            m_A_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_rawspecial(i);
                })
            &&
            m_B_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G2, G1>& a) {
                    return a.marshal_in_rawspecial(i);
                })
            &&
            m_C_query.marshal_in(
                is,
                [] (std::istream& i, Pairing<G1, G1>& a) {
                    return a.marshal_in_rawspecial(i);
                })
            &&
            snarklib::marshal_in_rawspecial(is, m_H_query)
            &&
            snarklib::marshal_in_rawspecial(is, m_K_query);
    }

    void clear() {
        m_A_query.clear();
        m_B_query.clear();
        m_C_query.clear();
        m_H_query.clear();
        m_K_query.clear();
    }

    bool empty() const {
        return
            m_A_query.empty() ||
            m_B_query.empty() ||
            m_C_query.empty() ||
            m_H_query.empty() ||
            m_K_query.empty();
    }

private:
    SparseVector<Pairing<G1, G1>> m_A_query;
    SparseVector<Pairing<G2, G1>> m_B_query;
    SparseVector<Pairing<G1, G1>> m_C_query;
    std::vector<G1> m_H_query;
    std::vector<G1> m_K_query;
};

////////////////////////////////////////////////////////////////////////////////
// Verification key
//

template <typename PAIRING>
class PPZK_VerificationKey
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    PPZK_VerificationKey() = default;

    // copy semantics
    PPZK_VerificationKey(const G2& alphaA_g2,
                         const G1& alphaB_g1,
                         const G2& alphaC_g2,
                         const G2& gamma_g2,
                         const G1& gamma_beta_g1,
                         const G2& gamma_beta_g2,
                         const G2& rC_Z_g2,
                         const PPZK_QueryIC<PAIRING>& encoded_IC_query)
        : m_alphaA_g2(alphaA_g2),
          m_alphaB_g1(alphaB_g1),
          m_alphaC_g2(alphaC_g2),
          m_gamma_g2(gamma_g2),
          m_gamma_beta_g1(gamma_beta_g1),
          m_gamma_beta_g2(gamma_beta_g2),
          m_rC_Z_g2(rC_Z_g2),
          m_encoded_IC_query(encoded_IC_query)
    {
        m_alphaA_g2.toSpecial();
        m_alphaB_g1.toSpecial();
        m_alphaC_g2.toSpecial();
        m_gamma_g2.toSpecial();
        m_gamma_beta_g1.toSpecial();
        m_gamma_beta_g2.toSpecial();
        m_rC_Z_g2.toSpecial();
        m_encoded_IC_query.toSpecial();
    }

    // move semantics
    PPZK_VerificationKey(const G2& alphaA_g2,
                         const G1& alphaB_g1,
                         const G2& alphaC_g2,
                         const G2& gamma_g2,
                         const G1& gamma_beta_g1,
                         const G2& gamma_beta_g2,
                         const G2& rC_Z_g2,
                         PPZK_QueryIC<PAIRING>&& encoded_IC_query)
        : m_alphaA_g2(alphaA_g2),
          m_alphaB_g1(alphaB_g1),
          m_alphaC_g2(alphaC_g2),
          m_gamma_g2(gamma_g2),
          m_gamma_beta_g1(gamma_beta_g1),
          m_gamma_beta_g2(gamma_beta_g2),
          m_rC_Z_g2(rC_Z_g2),
          m_encoded_IC_query(std::move(encoded_IC_query))
    {
        m_alphaA_g2.toSpecial();
        m_alphaB_g1.toSpecial();
        m_alphaC_g2.toSpecial();
        m_gamma_g2.toSpecial();
        m_gamma_beta_g1.toSpecial();
        m_gamma_beta_g2.toSpecial();
        m_rC_Z_g2.toSpecial();
        m_encoded_IC_query.toSpecial();
    }

    const G2& alphaA_g2() const { return m_alphaA_g2; }
    const G1& alphaB_g1() const { return m_alphaB_g1; }
    const G2& alphaC_g2() const { return m_alphaC_g2; }
    const G2& gamma_g2() const { return m_gamma_g2; }
    const G1& gamma_beta_g1() const { return m_gamma_beta_g1; }
    const G2& gamma_beta_g2() const { return m_gamma_beta_g2; }
    const G2& rC_Z_g2() const { return m_rC_Z_g2; }

    const PPZK_QueryIC<PAIRING>& encoded_IC_query() const {
        return m_encoded_IC_query;
    }

    bool operator== (const PPZK_VerificationKey& other) const {
        return
            alphaA_g2() == other.alphaA_g2() &&
            alphaB_g1() == other.alphaB_g1() &&
            alphaC_g2() == other.alphaC_g2() &&
            gamma_g2() == other.gamma_g2() &&
            gamma_beta_g1() == other.gamma_beta_g1() &&
            gamma_beta_g2() == other.gamma_beta_g2() &&
            rC_Z_g2() == other.rC_Z_g2() &&
            encoded_IC_query() == other.encoded_IC_query();
    }

    bool operator!= (const PPZK_VerificationKey& other) const {
        return ! (*this == other);
    }

    void marshal_out(std::ostream& os) const {
        alphaA_g2().marshal_out(os);
        alphaB_g1().marshal_out(os);
        alphaC_g2().marshal_out(os);
        gamma_g2().marshal_out(os);
        gamma_beta_g1().marshal_out(os);
        gamma_beta_g2().marshal_out(os);
        rC_Z_g2().marshal_out(os);
        encoded_IC_query().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_alphaA_g2.marshal_in(is) &&
            m_alphaB_g1.marshal_in(is) &&
            m_alphaC_g2.marshal_in(is) &&
            m_gamma_g2.marshal_in(is) &&
            m_gamma_beta_g1.marshal_in(is) &&
            m_gamma_beta_g2.marshal_in(is) &&
            m_rC_Z_g2.marshal_in(is) &&
            m_encoded_IC_query.marshal_in(is);
    }

    void marshal_out_raw(std::ostream& os) const {
        alphaA_g2().marshal_out_raw(os);
        alphaB_g1().marshal_out_raw(os);
        alphaC_g2().marshal_out_raw(os);
        gamma_g2().marshal_out_raw(os);
        gamma_beta_g1().marshal_out_raw(os);
        gamma_beta_g2().marshal_out_raw(os);
        rC_Z_g2().marshal_out_raw(os);
        encoded_IC_query().marshal_out_raw(os);
    }

    bool marshal_in_raw(std::istream& is) {
        return
            m_alphaA_g2.marshal_in_raw(is) &&
            m_alphaB_g1.marshal_in_raw(is) &&
            m_alphaC_g2.marshal_in_raw(is) &&
            m_gamma_g2.marshal_in_raw(is) &&
            m_gamma_beta_g1.marshal_in_raw(is) &&
            m_gamma_beta_g2.marshal_in_raw(is) &&
            m_rC_Z_g2.marshal_in_raw(is) &&
            m_encoded_IC_query.marshal_in_raw(is);
    }

    void marshal_out_special(std::ostream& os) const {
        alphaA_g2().marshal_out_special(os);
        alphaB_g1().marshal_out_special(os);
        alphaC_g2().marshal_out_special(os);
        gamma_g2().marshal_out_special(os);
        gamma_beta_g1().marshal_out_special(os);
        gamma_beta_g2().marshal_out_special(os);
        rC_Z_g2().marshal_out_special(os);
        encoded_IC_query().marshal_out_special(os);
    }

    bool marshal_in_special(std::istream& is) {
        return
            m_alphaA_g2.marshal_in_special(is) &&
            m_alphaB_g1.marshal_in_special(is) &&
            m_alphaC_g2.marshal_in_special(is) &&
            m_gamma_g2.marshal_in_special(is) &&
            m_gamma_beta_g1.marshal_in_special(is) &&
            m_gamma_beta_g2.marshal_in_special(is) &&
            m_rC_Z_g2.marshal_in_special(is) &&
            m_encoded_IC_query.marshal_in_special(is);
    }

    void marshal_out_rawspecial(std::ostream& os) const {
        alphaA_g2().marshal_out_rawspecial(os);
        alphaB_g1().marshal_out_rawspecial(os);
        alphaC_g2().marshal_out_rawspecial(os);
        gamma_g2().marshal_out_rawspecial(os);
        gamma_beta_g1().marshal_out_rawspecial(os);
        gamma_beta_g2().marshal_out_rawspecial(os);
        rC_Z_g2().marshal_out_rawspecial(os);
        encoded_IC_query().marshal_out_rawspecial(os);
    }

    bool marshal_in_rawspecial(std::istream& is) {
        return
            m_alphaA_g2.marshal_in_rawspecial(is) &&
            m_alphaB_g1.marshal_in_rawspecial(is) &&
            m_alphaC_g2.marshal_in_rawspecial(is) &&
            m_gamma_g2.marshal_in_rawspecial(is) &&
            m_gamma_beta_g1.marshal_in_rawspecial(is) &&
            m_gamma_beta_g2.marshal_in_rawspecial(is) &&
            m_rC_Z_g2.marshal_in_rawspecial(is) &&
            m_encoded_IC_query.marshal_in_rawspecial(is);
    }

    void clear() {
        m_alphaA_g2 = G2::zero();
        m_alphaB_g1 = G1::zero();
        m_alphaC_g2 = G2::zero();
        m_gamma_g2 = G2::zero();
        m_gamma_beta_g1 = G1::zero();
        m_gamma_beta_g2 = G2::zero();
        m_rC_Z_g2 = G2::zero();
        m_encoded_IC_query.clear();
    }

    bool empty() const {
        return
            m_alphaA_g2.isZero() ||
            m_alphaB_g1.isZero() ||
            m_alphaC_g2.isZero() ||
            m_gamma_g2.isZero() ||
            m_gamma_beta_g1.isZero() ||
            m_gamma_beta_g2.isZero() ||
            m_rC_Z_g2.isZero() ||
            m_encoded_IC_query.empty();
    }

private:
    G2 m_alphaA_g2;
    G1 m_alphaB_g1;
    G2 m_alphaC_g2;
    G2 m_gamma_g2;
    G1 m_gamma_beta_g1;
    G2 m_gamma_beta_g2;
    G2 m_rC_Z_g2;
    PPZK_QueryIC<PAIRING> m_encoded_IC_query;
};

////////////////////////////////////////////////////////////////////////////////
// Precomputed verification key (Miller loop input)
//

template <typename PAIRING>
class PPZK_PrecompVerificationKey
{
    typedef typename PAIRING::G2 G2;
    typedef typename PAIRING::G1_precomp G1_precomp;
    typedef typename PAIRING::G2_precomp G2_precomp;

public:
    PPZK_PrecompVerificationKey(const PPZK_VerificationKey<PAIRING>& vk)
        : m_pp_G2_one_precomp(G2::one()),
          m_vk_alphaA_g2_precomp(vk.alphaA_g2()),
          m_vk_alphaB_g1_precomp(vk.alphaB_g1()),
          m_vk_alphaC_g2_precomp(vk.alphaC_g2()),
          m_vk_rC_Z_g2_precomp(vk.rC_Z_g2()),
          m_vk_gamma_g2_precomp(vk.gamma_g2()),
          m_vk_gamma_beta_g1_precomp(vk.gamma_beta_g1()),
          m_vk_gamma_beta_g2_precomp(vk.gamma_beta_g2()),
          m_encoded_IC_query(vk.encoded_IC_query())
    {}

    const G2_precomp& pp_G2_one_precomp() const { return m_pp_G2_one_precomp; }
    const G2_precomp& vk_alphaA_g2_precomp() const { return m_vk_alphaA_g2_precomp; }
    const G1_precomp& vk_alphaB_g1_precomp() const { return m_vk_alphaB_g1_precomp; }
    const G2_precomp& vk_alphaC_g2_precomp() const { return m_vk_alphaC_g2_precomp; }
    const G2_precomp& vk_rC_Z_g2_precomp() const { return m_vk_rC_Z_g2_precomp; }
    const G2_precomp& vk_gamma_g2_precomp() const { return m_vk_gamma_g2_precomp; }
    const G1_precomp& vk_gamma_beta_g1_precomp() const { return m_vk_gamma_beta_g1_precomp; }
    const G2_precomp& vk_gamma_beta_g2_precomp() const { return m_vk_gamma_beta_g2_precomp; }

    const PPZK_QueryIC<PAIRING>& encoded_IC_query() const {
        return m_encoded_IC_query;
    }

private:
    G2_precomp m_pp_G2_one_precomp;
    G2_precomp m_vk_alphaA_g2_precomp;
    G1_precomp m_vk_alphaB_g1_precomp;
    G2_precomp m_vk_alphaC_g2_precomp;
    G2_precomp m_vk_rC_Z_g2_precomp;
    G2_precomp m_vk_gamma_g2_precomp;
    G1_precomp m_vk_gamma_beta_g1_precomp;
    G2_precomp m_vk_gamma_beta_g2_precomp;
    PPZK_QueryIC<PAIRING> m_encoded_IC_query;
};

} // namespace snarklib

#endif
