#ifndef _SNARKLIB_AUTOTEST_PPZK_HPP_
#define _SNARKLIB_AUTOTEST_PPZK_HPP_

#include <cstdint>
#include <vector>
#include "AutoTest.hpp"
#include "AutoTest_R1CS.hpp"
#include "AuxSTL.hpp"
#include "common/types.hpp"
#include "Pairing.hpp"
#include "PPZK.hpp"
#include "r1cs_ppzksnark/r1cs_ppzksnark.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// end-to-end using libsnark
//

template <typename T, typename U>
class AutoTest_PPZK_libsnark_only : public AutoTest
{
    typedef typename libsnark::default_pp PPT;

public:
    AutoTest_PPZK_libsnark_only(const AutoTestR1CS<T, U>& cs)
        : AutoTest(cs),
          m_cs(cs.systemA()),
          m_witness(cs.witnessA()),
          m_input(cs.inputA())
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_cs);

        const auto pvk
            = libsnark::r1cs_ppzksnark_verifier_process_vk<PPT>(
                keypair.vk);

        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_witness);

        const auto ans_online
            = libsnark::r1cs_ppzksnark_online_verifier_strong_IC<PPT>(
                pvk,
                m_input,
                proof);

        checkPass(ans_online);
    }

private:
    const libsnark::r1cs_constraint_system<U> m_cs;
    const libsnark::r1cs_variable_assignment<U> m_witness;
    const libsnark::r1cs_variable_assignment<U> m_input;
};

////////////////////////////////////////////////////////////////////////////////
// verification uses redesigned code
//

template <typename PAIRING, typename U>
class AutoTest_PPZK_strongVerify : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    typedef typename libsnark::default_pp PPT;

public:
    AutoTest_PPZK_strongVerify(const AutoTestR1CS<Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_constraintSystem.systemA());

        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.witnessA());

        // encoded IC query
        G1 base;
        const size_t encSize = keypair.vk.encoded_IC_query->encoded_terms.size();
        std::vector<G1> encoded_terms(encSize);
        copyData(keypair.vk.encoded_IC_query->base, base);
        for (std::size_t i = 0; i < encSize; ++i) {
            copyData(keypair.vk.encoded_IC_query->encoded_terms[i], encoded_terms[i]);
        }
        const PPZK_IC_Query<PAIRING> icqB(base, encoded_terms);

        // verification key
        G1 alphaB_g1, gamma_beta_g1;
        G2 alphaA_g2, alphaC_g2, gamma_g2, gamma_beta_g2, rC_Z_g2;
        copyData(keypair.vk.alphaA_g2, alphaA_g2);
        copyData(keypair.vk.alphaB_g1, alphaB_g1);
        copyData(keypair.vk.alphaC_g2, alphaC_g2);
        copyData(keypair.vk.gamma_g2, gamma_g2);
        copyData(keypair.vk.gamma_beta_g1, gamma_beta_g1);
        copyData(keypair.vk.gamma_beta_g2, gamma_beta_g2);
        copyData(keypair.vk.rC_Z_g2, rC_Z_g2);
        const PPZK_VerificationKey<PAIRING> vkB(alphaA_g2,
                                                alphaB_g1,
                                                alphaC_g2,
                                                gamma_g2,
                                                gamma_beta_g1,
                                                gamma_beta_g2,
                                                rC_Z_g2,
                                                icqB);

        // proof
        G1 AG, AH, BH, CG, CH, H, K;
        G2 BG;
        copyData(proof.g_A.g, AG);
        copyData(proof.g_A.h, AH);
        copyData(proof.g_B.g, BG);
        copyData(proof.g_B.h, BH);
        copyData(proof.g_C.g, CG);
        copyData(proof.g_C.h, CH);
        copyData(proof.g_H, H);
        copyData(proof.g_K, K);
        const PPZK_Proof<PAIRING> proofB(Pairing<G1, G1>(AG, AH),
                                         Pairing<G2, G1>(BG, BH),
                                         Pairing<G1, G1>(CG, CH),
                                         H,
                                         K);

        const auto ans
            = strongVerify(
                vkB,
                m_constraintSystem.inputB(),
                proofB);

        checkPass(ans);
    }

private:
    const AutoTestR1CS<Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// compare original and redesigned proof output
//

template <typename PAIRING, typename U>
class AutoTest_PPZK_ProofCompare : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    typedef typename libsnark::default_pp PPT;

public:
    AutoTest_PPZK_ProofCompare(const AutoTestR1CS<Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_constraintSystem.systemA());

        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.witnessA());

        // proof from original code
        G1 AG, AH, BH, CG, CH, H, K;
        G2 BG;
        copyData(proof.g_A.g, AG);
        copyData(proof.g_A.h, AH);
        copyData(proof.g_B.g, BG);
        copyData(proof.g_B.h, BH);
        copyData(proof.g_C.g, CG);
        copyData(proof.g_C.h, CH);
        copyData(proof.g_H, H);
        copyData(proof.g_K, K);
        const PPZK_Proof<PAIRING> proofFromOriginal(
            Pairing<G1, G1>(AG, AH),
            Pairing<G2, G1>(BG, BH),
            Pairing<G1, G1>(CG, CH),
            H,
            K);

        // proving key
        SparseVector<Pairing<G1, G1>> A_query, C_query;
        SparseVector<Pairing<G2, G1>> B_query;
        std::vector<G1> H_query, K_query;
        copyData(keypair.pk.A_query, A_query);
        copyData(keypair.pk.B_query, B_query);
        copyData(keypair.pk.C_query, C_query);
        copyData(keypair.pk.H_query, H_query);
        copyData(keypair.pk.K_query, K_query);
        const PPZK_ProvingKey<PAIRING> pkB(A_query,
                                           B_query,
                                           C_query,
                                           H_query,
                                           K_query);

        // proof from redesigned code
        const PPZK_Proof<PAIRING> proofFromRedesign(
            m_constraintSystem.systemB(),
            m_constraintSystem.numberInputs(),
            pkB,
            m_constraintSystem.witnessB());

        // compare proofs (expect different results because of random numbers)
        if (! checkPass(proofFromOriginal.A() != proofFromRedesign.A())) {
            std::cout << "original A.G " << proofFromOriginal.A().G() << std::endl
                      << "redesign A.G " << proofFromRedesign.A().G() << std::endl
                      << "original A.H " << proofFromOriginal.A().H() << std::endl
                      << "redesign A.H " << proofFromRedesign.A().H() << std::endl;
        }

        if (! checkPass(proofFromOriginal.B() != proofFromRedesign.B())) {
            std::cout << "original B.G " << proofFromOriginal.B().G() << std::endl
                      << "redesign B.G " << proofFromRedesign.B().G() << std::endl
                      << "original B.H " << proofFromOriginal.B().H() << std::endl
                      << "redesign B.H " << proofFromRedesign.B().H() << std::endl;
        }

        if (! checkPass(proofFromOriginal.C() != proofFromRedesign.C())) {
            std::cout << "original C.G " << proofFromOriginal.C().G() << std::endl
                      << "redesign C.G " << proofFromRedesign.C().G() << std::endl
                      << "original C.H " << proofFromOriginal.C().H() << std::endl
                      << "redesign C.H " << proofFromRedesign.C().H() << std::endl;
        }

        if (! checkPass(proofFromOriginal.H() != proofFromRedesign.H())) {
            std::cout << "original H " << proofFromOriginal.H() << std::endl
                      << "redesign H " << proofFromRedesign.H() << std::endl;
        }

        if (! checkPass(proofFromOriginal.K() != proofFromRedesign.K())) {
            std::cout << "original K " << proofFromOriginal.K() << std::endl
                      << "redesign K " << proofFromRedesign.K() << std::endl;
        }
    }

private:
    const AutoTestR1CS<Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// verification and proof use redesigned code
//

template <typename PAIRING, typename U>
class AutoTest_PPZK_Proof : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    typedef typename libsnark::default_pp PPT;

public:
    AutoTest_PPZK_Proof(const AutoTestR1CS<Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_constraintSystem.systemA());

        // proving key
        SparseVector<Pairing<G1, G1>> A_query, C_query;
        SparseVector<Pairing<G2, G1>> B_query;
        std::vector<G1> H_query, K_query;
        copyData(keypair.pk.A_query, A_query);
        copyData(keypair.pk.B_query, B_query);
        copyData(keypair.pk.C_query, C_query);
        copyData(keypair.pk.H_query, H_query);
        copyData(keypair.pk.K_query, K_query);
        const PPZK_ProvingKey<PAIRING> pkB(A_query,
                                           B_query,
                                           C_query,
                                           H_query,
                                           K_query);

        // encoded IC query
        G1 base;
        const size_t encSize = keypair.vk.encoded_IC_query->encoded_terms.size();
        std::vector<G1> encoded_terms(encSize);
        copyData(keypair.vk.encoded_IC_query->base, base);
        for (std::size_t i = 0; i < encSize; ++i) {
            copyData(keypair.vk.encoded_IC_query->encoded_terms[i], encoded_terms[i]);
        }
        const PPZK_IC_Query<PAIRING> icqB(base, encoded_terms);

        // verification key
        G1 alphaB_g1, gamma_beta_g1;
        G2 alphaA_g2, alphaC_g2, gamma_g2, gamma_beta_g2, rC_Z_g2;
        copyData(keypair.vk.alphaA_g2, alphaA_g2);
        copyData(keypair.vk.alphaB_g1, alphaB_g1);
        copyData(keypair.vk.alphaC_g2, alphaC_g2);
        copyData(keypair.vk.gamma_g2, gamma_g2);
        copyData(keypair.vk.gamma_beta_g1, gamma_beta_g1);
        copyData(keypair.vk.gamma_beta_g2, gamma_beta_g2);
        copyData(keypair.vk.rC_Z_g2, rC_Z_g2);
        const PPZK_VerificationKey<PAIRING> vkB(alphaA_g2,
                                                alphaB_g1,
                                                alphaC_g2,
                                                gamma_g2,
                                                gamma_beta_g1,
                                                gamma_beta_g2,
                                                rC_Z_g2,
                                                icqB);

        // proof
        const PPZK_Proof<PAIRING> proofB(m_constraintSystem.systemB(),
                                         m_constraintSystem.numberInputs(),
                                         pkB,
                                         m_constraintSystem.witnessB());

        const auto ans
            = strongVerify(
                vkB,
                m_constraintSystem.inputB(),
                proofB);

        checkPass(ans);
    }

private:
    const AutoTestR1CS<Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// end-to-end using snarklib
//

template <typename PAIRING, typename U>
class AutoTest_PPZK_full_redesign : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    typedef typename libsnark::default_pp PPT;

public:
    AutoTest_PPZK_full_redesign(const AutoTestR1CS<Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const PPZK_Keypair<PAIRING> keypair(
            m_constraintSystem.systemB(),
            m_constraintSystem.numberInputs());

        const PPZK_Proof<PAIRING> proofB(m_constraintSystem.systemB(),
                                         m_constraintSystem.numberInputs(),
                                         keypair.pk(),
                                         m_constraintSystem.witnessB());

        const auto ans = strongVerify(keypair.vk(),
                                      m_constraintSystem.inputB(),
                                      proofB);

        checkPass(ans);
    }

private:
    const AutoTestR1CS<Fr, U> m_constraintSystem;
};

} // namespace snarklib

#endif
