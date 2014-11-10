#ifndef _SNARKLIB_PPZK_HPP_
#define _SNARKLIB_PPZK_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
#include "AuxSTL.hpp"
#include "Group.hpp"
#include "MultiExp.hpp"
#include "Pairing.hpp"
#include "QAP.hpp"
#include "Rank1DSL.hpp"
#include "WindowExp.hpp"

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
    {}

    const SparseVector<Pairing<G1, G1>>& A_query() const { return m_A_query; }
    const SparseVector<Pairing<G2, G1>>& B_query() const { return m_B_query; }
    const SparseVector<Pairing<G1, G1>>& C_query() const { return m_C_query; }
    const std::vector<G1>& H_query() const { return m_H_query; }
    const std::vector<G1>& K_query() const { return m_K_query; }

private:
    SparseVector<Pairing<G1, G1>> m_A_query;
    SparseVector<Pairing<G2, G1>> m_B_query;
    SparseVector<Pairing<G1, G1>> m_C_query;
    std::vector<G1> m_H_query;
    std::vector<G1> m_K_query;
};

////////////////////////////////////////////////////////////////////////////////
// Input consistency
//

template <typename PAIRING>
class PPZK_IC_Query
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;

public:
    PPZK_IC_Query() = default;

    PPZK_IC_Query(const G1& base,
                  const std::vector<G1>& encoded_terms)
        : m_base(base),
          m_encoded_terms(encoded_terms)
    {}

    PPZK_IC_Query(const WindowExp<G1>& g1_table,
                  const std::vector<Fr>& coeffs)
        : PPZK_IC_Query{
            coeffs[0] * G1::one(),
            g1_table.batchExp(
                std::vector<Fr>(coeffs.begin() + 1, coeffs.end()))}
    {}

    PPZK_IC_Query accumulate(const R1Witness<Fr>& witness) const
    {
        G1 base = m_base;
        std::vector<G1> encoded_terms;

        const std::size_t
            wsize = witness.size(),
            tsize = m_encoded_terms.size();

        if (wsize < tsize) {
            base = base + multiExp(
                std::vector<G1>(m_encoded_terms.begin(),
                                m_encoded_terms.begin() + wsize),
                *witness);

            encoded_terms
                = std::vector<G1>(m_encoded_terms.begin() + wsize,
                                  m_encoded_terms.end());

        } else if (wsize > tsize) {
            base = base + multiExp(m_encoded_terms,
                                   *witness.truncate(tsize));

        } else {
            base = base + multiExp(m_encoded_terms,
                                   *witness);
        }

        return PPZK_IC_Query(base, encoded_terms);
    }

    const G1& base() const {
        return m_base;
    }

    std::size_t input_size() const {
        return m_encoded_terms.size();
    }

private:
    G1 m_base;
    std::vector<G1> m_encoded_terms;
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

    PPZK_VerificationKey(const G2& alphaA_g2,
                         const G1& alphaB_g1,
                         const G2& alphaC_g2,
                         const G2& gamma_g2,
                         const G1& gamma_beta_g1,
                         const G2& gamma_beta_g2,
                         const G2& rC_Z_g2,
                         const PPZK_IC_Query<PAIRING>& encoded_IC_query)
        : m_alphaA_g2(alphaA_g2),
          m_alphaB_g1(alphaB_g1),
          m_alphaC_g2(alphaC_g2),
          m_gamma_g2(gamma_g2),
          m_gamma_beta_g1(gamma_beta_g1),
          m_gamma_beta_g2(gamma_beta_g2),
          m_rC_Z_g2(rC_Z_g2),
          m_encoded_IC_query(encoded_IC_query)
    {}

    const G2& alphaA_g2() const { return m_alphaA_g2; }
    const G1& alphaB_g1() const { return m_alphaB_g1; }
    const G2& alphaC_g2() const { return m_alphaC_g2; }
    const G2& gamma_g2() const { return m_gamma_g2; }
    const G1& gamma_beta_g1() const { return m_gamma_beta_g1; }
    const G2& gamma_beta_g2() const { return m_gamma_beta_g2; }
    const G2& rC_Z_g2() const { return m_rC_Z_g2; }

    const PPZK_IC_Query<PAIRING>& encoded_IC_query() const {
        return m_encoded_IC_query;
    }

private:
    G2 m_alphaA_g2;
    G1 m_alphaB_g1;
    G2 m_alphaC_g2;
    G2 m_gamma_g2;
    G1 m_gamma_beta_g1;
    G2 m_gamma_beta_g2;
    G2 m_rC_Z_g2;
    PPZK_IC_Query<PAIRING> m_encoded_IC_query;
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

    const PPZK_IC_Query<PAIRING>& encoded_IC_query() const {
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
    PPZK_IC_Query<PAIRING> m_encoded_IC_query;
};

////////////////////////////////////////////////////////////////////////////////
// Key pair: proving and verification
//

template <typename PAIRING>
class PPZK_Keypair
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    PPZK_Keypair(const R1System<Fr>& constraintSystem,
                 const std::size_t numCircuitInputs)
    {
        const auto
            point = Fr::random(),
            alphaA = Fr::random(),
            alphaB = Fr::random(),
            alphaC = Fr::random(),
            rA = Fr::random(),
            rB = Fr::random(),
            beta = Fr::random(),
            gamma = Fr::random();

        const auto rC = rA * rB;

        QAP_ABCH<Fr> abch(constraintSystem,
                          numCircuitInputs,
                          point);

        const auto Z = abch.FFT()->compute_Z(point);

        const WindowExp<G1> g1_table(WindowExp<G1>::windowBits(abch.g1_exp_count()));
        const WindowExp<G2> g2_table(WindowExp<G2>::windowBits(abch.g2_exp_count()));

        auto K_query = g1_table.batchExp(abch.K_query(beta, rA, rB, rC));
#ifdef USE_ADD_SPECIAL
        batchSpecial(K_query);
#endif

        // zero out IC from A_query and place it into IC coefficients
        const auto IC_coefficients = abch.IC_coefficients(rA);

        m_pk = PPZK_ProvingKey<PAIRING>(
            batchExp(g1_table, g1_table, rA, rA * alphaA, abch.A_query()),
            batchExp(g2_table, g1_table, rB, rB * alphaB, abch.B_query()),
            batchExp(g1_table, g1_table, rC, rC * alphaC, abch.C_query()),
            g1_table.batchExp(abch.H_query()),
            K_query);

        m_vk = PPZK_VerificationKey<PAIRING>(
            alphaA * G2::one(),
            alphaB * G1::one(),
            alphaC * G2::one(),
            gamma * G2::one(),
            (gamma * beta) * G1::one(),
            (gamma * beta) * G2::one(),
            (rC * Z) * G2::one(),
            PPZK_IC_Query<PAIRING>(g1_table, IC_coefficients));
    }

    const PPZK_ProvingKey<PAIRING>& pk() const { return m_pk; }
    const PPZK_VerificationKey<PAIRING>& vk() const { return m_vk; }

private:
    PPZK_ProvingKey<PAIRING> m_pk;
    PPZK_VerificationKey<PAIRING> m_vk;
};

////////////////////////////////////////////////////////////////////////////////
// Proof
//

template <typename PAIRING>
class PPZK_Proof
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    PPZK_Proof(const Pairing<G1, G1>& A,
               const Pairing<G2, G1>& B,
               const Pairing<G1, G1>& C,
               const G1& H,
               const G1& K)
        : m_A(A),
          m_B(B),
          m_C(C),
          m_H(H),
          m_K(K)
    {}

    PPZK_Proof(const R1System<Fr>& constraintSystem,
               const std::size_t numCircuitInputs,
               const PPZK_ProvingKey<PAIRING>& pk,
               const R1Witness<Fr>& witness)
    {
        const auto
            d1 = Fr::random(),
            d2 = Fr::random(),
            d3 = Fr::random();

        const QAP_Witness<Fr> qap(
            constraintSystem, numCircuitInputs, witness, d1, d2, d3);

        const auto& A_query = pk.A_query();
        const auto& B_query = pk.B_query();
        const auto& C_query = pk.C_query();
        const auto& H_query = pk.H_query();
        const auto& K_query = pk.K_query();

        m_A = (d1 * A_query.getElementForIndex(0))
            + A_query.getElementForIndex(3)
            + multiExp01(A_query, *witness, 4, 4 + qap.numVariables());

        m_B = (d2 * B_query.getElementForIndex(1))
            + B_query.getElementForIndex(3)
            + multiExp01(B_query, *witness, 4, 4 + qap.numVariables());

        m_C = (d3 * C_query.getElementForIndex(2))
            + C_query.getElementForIndex(3)
            + multiExp01(C_query, *witness, 4, 4 + qap.numVariables());

        m_H = multiExp(H_query, qap.H());

        m_K = (d1 * K_query[0]) + (d2 * K_query[1]) + (d3 * K_query[2])
            + K_query[3]
            + multiExp01(
                std::vector<G1>(K_query.begin() + 4, K_query.end()),
                *witness);
    }

    const Pairing<G1, G1>& A() const { return m_A; }
    const Pairing<G2, G1>& B() const { return m_B; }
    const Pairing<G1, G1>& C() const { return m_C; }
    const G1& H() const { return m_H; }
    const G1& K() const { return m_K; }

    bool wellFormed() const {
        return
            m_A.G().wellFormed() && m_A.H().wellFormed() &&
            m_B.G().wellFormed() && m_B.H().wellFormed() &&
            m_C.G().wellFormed() && m_C.H().wellFormed() &&
            m_H.wellFormed() &&
            m_K.wellFormed();
    }

private:
    Pairing<G1, G1> m_A;
    Pairing<G2, G1> m_B;
    Pairing<G1, G1> m_C;
    G1 m_H;
    G1 m_K;
};

////////////////////////////////////////////////////////////////////////////////
// Verification functions
//

template <typename PAIRING>
bool weakVerify(const PPZK_PrecompVerificationKey<PAIRING>& pvk,
                const R1Witness<typename PAIRING::Fr>& input,
                const PPZK_Proof<PAIRING>& proof)
{
    typedef typename PAIRING::GT GT;
    typedef typename PAIRING::G1_precomp G1_precomp;
    typedef typename PAIRING::G2_precomp G2_precomp;

    const auto accum_IC = pvk.encoded_IC_query().accumulate(input);
    assert(0 == accum_IC.input_size());

    if (! proof.wellFormed()) {
        return false;
    }

    const auto ONE = GT::one();

     // knowledge commitment for A
    const G1_precomp proof_g_A_g_precomp(proof.A().G());
    const G1_precomp proof_g_A_h_precomp(proof.A().H());
    const auto kc_A_1 = PAIRING::ate_miller_loop(proof_g_A_g_precomp, pvk.vk_alphaA_g2_precomp());
    const auto kc_A_2 = PAIRING::ate_miller_loop(proof_g_A_h_precomp, pvk.pp_G2_one_precomp());
    const auto kc_A = PAIRING::final_exponentiation(kc_A_1 * unitary_inverse(kc_A_2));
    if (ONE != kc_A) {
        return false;
    }
    
    // knowledge commitment for B
    const G2_precomp proof_g_B_g_precomp(proof.B().G());
    const G1_precomp proof_g_B_h_precomp(proof.B().H());
    const auto kc_B_1 = PAIRING::ate_miller_loop(pvk.vk_alphaB_g1_precomp(), proof_g_B_g_precomp);
    const auto kc_B_2 = PAIRING::ate_miller_loop(proof_g_B_h_precomp, pvk.pp_G2_one_precomp());
    const auto kc_B = PAIRING::final_exponentiation(kc_B_1 * unitary_inverse(kc_B_2));
    if (ONE != kc_B) {
        return false;
    }

    // knowledge commitment for C
    const G1_precomp proof_g_C_g_precomp(proof.C().G());
    const G1_precomp proof_g_C_h_precomp(proof.C().H());
    const auto kc_C_1 = PAIRING::ate_miller_loop(proof_g_C_g_precomp, pvk.vk_alphaC_g2_precomp());
    const auto kc_C_2 = PAIRING::ate_miller_loop(proof_g_C_h_precomp, pvk.pp_G2_one_precomp());
    const auto kc_C = PAIRING::final_exponentiation(kc_C_1 * unitary_inverse(kc_C_2));
    if (ONE != kc_C) {
        return false;
    }

    // quadratic arithmetic program divisibility
    const G1_precomp proof_g_A_g_acc_precomp(proof.A().G() + accum_IC.base());
    const G1_precomp proof_g_H_precomp(proof.H());
    const auto QAP_1 = PAIRING::ate_miller_loop(proof_g_A_g_acc_precomp, proof_g_B_g_precomp);
    const auto QAP_23 = PAIRING::ate_double_miller_loop(proof_g_H_precomp, pvk.vk_rC_Z_g2_precomp(), proof_g_C_g_precomp, pvk.pp_G2_one_precomp());
    const auto QAP = PAIRING::final_exponentiation(QAP_1 * unitary_inverse(QAP_23));
    if (ONE != QAP) {
        return false;
    }

    // same coefficients
    const G1_precomp proof_g_K_precomp(proof.K());
    const G1_precomp proof_g_A_g_acc_C_precomp(proof.A().G() + accum_IC.base() + proof.C().G());
    const auto K_1 = PAIRING::ate_miller_loop(proof_g_K_precomp, pvk.vk_gamma_g2_precomp());
    const auto K_23 = PAIRING::ate_double_miller_loop(proof_g_A_g_acc_C_precomp, pvk.vk_gamma_beta_g2_precomp(), pvk.vk_gamma_beta_g1_precomp(), proof_g_B_g_precomp);
    const auto K = PAIRING::final_exponentiation(K_1 * unitary_inverse(K_23));
    if (ONE != K) {
        return false;
    }

    return true;
}

template <typename PAIRING>
bool weakVerify(const PPZK_VerificationKey<PAIRING>& vk,
                const R1Witness<typename PAIRING::Fr>& input,
                const PPZK_Proof<PAIRING>& proof)
{
    return weakVerify(PPZK_PrecompVerificationKey<PAIRING>(vk),
                      input,
                      proof);
}

template <typename PAIRING>
bool strongVerify(const PPZK_PrecompVerificationKey<PAIRING>& pvk,
                  const R1Witness<typename PAIRING::Fr>& input,
                  const PPZK_Proof<PAIRING>& proof)
{
    return (pvk.encoded_IC_query().input_size() == input.size())
        ? weakVerify(pvk, input, proof)
        : false;
}

template <typename PAIRING>
bool strongVerify(const PPZK_VerificationKey<PAIRING>& vk,
                  const R1Witness<typename PAIRING::Fr>& input,
                  const PPZK_Proof<PAIRING>& proof)
{
    return strongVerify(PPZK_PrecompVerificationKey<PAIRING>(vk),
                        input,
                        proof);
}


} // namespace snarklib

#endif
