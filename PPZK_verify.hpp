#ifndef _SNARKLIB_PPZK_VERIFY_HPP_
#define _SNARKLIB_PPZK_VERIFY_HPP_

#include <memory>

#include <snarklib/ProgressCallback.hpp>
#include <snarklib/PPZK_keystruct.hpp>
#include <snarklib/PPZK_proof.hpp>
#include <snarklib/Rank1DSL.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Verification functions
//

template <typename PAIRING>
bool weakVerify(const PPZK_PrecompVerificationKey<PAIRING>& pvk,
                const R1Witness<typename PAIRING::Fr>& input,
                const PPZK_Proof<PAIRING>& proof,
                ProgressCallback* callback = nullptr)
{
    ProgressCallback_NOP<PAIRING> dummyNOP;
    ProgressCallback* dummy = callback ? callback : std::addressof(dummyNOP);
    dummy->majorSteps(6);

    typedef typename PAIRING::GT GT;
    typedef typename PAIRING::G1_precomp G1_precomp;
    typedef typename PAIRING::G2_precomp G2_precomp;

    const auto ONE = GT::one();

    // step 6 (starting) - accumulate input consistency
    dummy->major();
    const auto accum_IC = pvk.encoded_IC_query().accumWitness(input);
    if (0 != accum_IC.input_size() || ! proof.wellFormed()) return false;

    // step 5 - knowledge commitment for A
    dummy->major();
    const G1_precomp proof_g_A_g_precomp(proof.A().G());
    const G1_precomp proof_g_A_h_precomp(proof.A().H());

    const auto kc_A_1 = PAIRING::ate_miller_loop(
        proof_g_A_g_precomp,
        pvk.vk_alphaA_g2_precomp());

    const auto kc_A_2 = PAIRING::ate_miller_loop(
        proof_g_A_h_precomp,
        pvk.pp_G2_one_precomp());

    const auto kc_A = PAIRING::final_exponentiation(
        kc_A_1 * unitary_inverse(kc_A_2));

    if (ONE != kc_A) return false;

    // step 4 - knowledge commitment for B
    dummy->major();
    const G2_precomp proof_g_B_g_precomp(proof.B().G());
    const G1_precomp proof_g_B_h_precomp(proof.B().H());

    const auto kc_B_1 = PAIRING::ate_miller_loop(
        pvk.vk_alphaB_g1_precomp(),
        proof_g_B_g_precomp);

    const auto kc_B_2 = PAIRING::ate_miller_loop(
        proof_g_B_h_precomp,
        pvk.pp_G2_one_precomp());

    const auto kc_B = PAIRING::final_exponentiation(
        kc_B_1 * unitary_inverse(kc_B_2));

    if (ONE != kc_B) return false;

    // step 3 - knowledge commitment for C
    dummy->major();
    const G1_precomp proof_g_C_g_precomp(proof.C().G());
    const G1_precomp proof_g_C_h_precomp(proof.C().H());

    const auto kc_C_1 = PAIRING::ate_miller_loop(
        proof_g_C_g_precomp,
        pvk.vk_alphaC_g2_precomp());

    const auto kc_C_2 = PAIRING::ate_miller_loop(
        proof_g_C_h_precomp,
        pvk.pp_G2_one_precomp());

    const auto kc_C = PAIRING::final_exponentiation(
        kc_C_1 * unitary_inverse(kc_C_2));

    if (ONE != kc_C) return false;

    // step 2 - quadratic arithmetic program divisibility
    dummy->major();
    const G1_precomp proof_g_A_g_acc_precomp(proof.A().G() + accum_IC.base());
    const G1_precomp proof_g_H_precomp(proof.H());

    const auto QAP_1 = PAIRING::ate_miller_loop(
        proof_g_A_g_acc_precomp,
        proof_g_B_g_precomp);

    const auto QAP_23 = PAIRING::ate_double_miller_loop(
        proof_g_H_precomp,
        pvk.vk_rC_Z_g2_precomp(),
        proof_g_C_g_precomp,
        pvk.pp_G2_one_precomp());

    const auto QAP = PAIRING::final_exponentiation(
        QAP_1 * unitary_inverse(QAP_23));

    if (ONE != QAP) return false;

    // step 1 - same coefficients
    dummy->major();
    const G1_precomp proof_g_K_precomp(proof.K());
    const G1_precomp proof_g_A_g_acc_C_precomp(proof.A().G() + accum_IC.base() + proof.C().G());

    const auto K_1 = PAIRING::ate_miller_loop(
        proof_g_K_precomp,
        pvk.vk_gamma_g2_precomp());

    const auto K_23 = PAIRING::ate_double_miller_loop(
        proof_g_A_g_acc_C_precomp,
        pvk.vk_gamma_beta_g2_precomp(),
        pvk.vk_gamma_beta_g1_precomp(),
        proof_g_B_g_precomp);

    const auto K = PAIRING::final_exponentiation(
        K_1 * unitary_inverse(K_23));

    if (ONE != K) return false;

    return true;
}

template <typename PAIRING>
bool weakVerify(const PPZK_VerificationKey<PAIRING>& vk,
                const R1Witness<typename PAIRING::Fr>& input,
                const PPZK_Proof<PAIRING>& proof,
                ProgressCallback* callback = nullptr)
{
    return weakVerify(PPZK_PrecompVerificationKey<PAIRING>(vk),
                      input,
                      proof,
                      callback);
}

template <typename PAIRING>
bool strongVerify(const PPZK_PrecompVerificationKey<PAIRING>& pvk,
                  const R1Witness<typename PAIRING::Fr>& input,
                  const PPZK_Proof<PAIRING>& proof,
                  ProgressCallback* callback = nullptr)
{
    return (pvk.encoded_IC_query().input_size() == input.size())
        ? weakVerify(pvk, input, proof, callback)
        : false;
}

template <typename PAIRING>
bool strongVerify(const PPZK_VerificationKey<PAIRING>& vk,
                  const R1Witness<typename PAIRING::Fr>& input,
                  const PPZK_Proof<PAIRING>& proof,
                  ProgressCallback* callback = nullptr)
{
    return strongVerify(PPZK_PrecompVerificationKey<PAIRING>(vk),
                        input,
                        proof,
                        callback);
}

} // namespace snarklib

#endif
