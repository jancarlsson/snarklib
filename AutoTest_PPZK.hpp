#ifndef _SNARKLIB_AUTOTEST_PPZK_HPP_
#define _SNARKLIB_AUTOTEST_PPZK_HPP_

#include <cstdint>
#include <vector>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "common/types.hpp"
#include /*libsnark*/ "r1cs_ppzksnark/r1cs_ppzksnark.hpp"
#else
#include /*libsnark*/ "algebra/curves/public_params.hpp"
#include /*libsnark*/ "zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/AutoTest_R1CS.hpp"
#include "snarklib/AuxSTL.hpp"
#include "snarklib/ForeignLib.hpp"
#include "snarklib/Pairing.hpp"
#include "snarklib/PPZK_keypair.hpp"
#include "snarklib/PPZK_keystruct.hpp"
#include "snarklib/PPZK_query.hpp"
#include "snarklib/PPZK_proof.hpp"
#include "snarklib/PPZK_verify.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// end-to-end using libsnark
//

template <typename T, typename U>
class AutoTest_PPZK_libsnark_only : public AutoTest
{
#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp PPT;
#else
    typedef libsnark::default_ec_pp PPT;
#endif

public:
    template <template <typename> class SYS>
    AutoTest_PPZK_libsnark_only(const AutoTestR1CS<SYS, T, U>& cs)
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

#ifdef USE_OLD_LIBSNARK
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_witness);
#else
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_input,
                m_witness);
#endif

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

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_strongVerify : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp PPT;
#else
    typedef libsnark::default_ec_pp PPT;
#endif

public:
    AutoTest_PPZK_strongVerify(const AutoTestR1CS<SYS, Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_constraintSystem.systemA());

#ifdef USE_OLD_LIBSNARK
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.witnessA());
#else
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.inputA(),
                m_constraintSystem.witnessA());
#endif

        // encoded IC query
        G1 base;
#ifdef USE_OLD_LIBSNARK
        const size_t encSize = keypair.vk.encoded_IC_query->encoded_terms.size();
#else
        const size_t encSize = keypair.vk.encoded_IC_query.rest.size();
#endif
        std::vector<G1> encoded_terms(encSize);
#ifdef USE_OLD_LIBSNARK
        copy_libsnark(keypair.vk.encoded_IC_query->base, base);
#else
        copy_libsnark(keypair.vk.encoded_IC_query.first, base);
#endif
        for (std::size_t i = 0; i < encSize; ++i) {
#ifdef USE_OLD_LIBSNARK
            copy_libsnark(keypair.vk.encoded_IC_query->encoded_terms[i], encoded_terms[i]);
#else
            copy_libsnark(keypair.vk.encoded_IC_query.rest[i], encoded_terms[i]);
#endif
        }
        const PPZK_QueryIC<PAIRING> icqB(base, encoded_terms);

        // verification key
        G1 alphaB_g1, gamma_beta_g1;
        G2 alphaA_g2, alphaC_g2, gamma_g2, gamma_beta_g2, rC_Z_g2;
        copy_libsnark(keypair.vk.alphaA_g2, alphaA_g2);
        copy_libsnark(keypair.vk.alphaB_g1, alphaB_g1);
        copy_libsnark(keypair.vk.alphaC_g2, alphaC_g2);
        copy_libsnark(keypair.vk.gamma_g2, gamma_g2);
        copy_libsnark(keypair.vk.gamma_beta_g1, gamma_beta_g1);
        copy_libsnark(keypair.vk.gamma_beta_g2, gamma_beta_g2);
        copy_libsnark(keypair.vk.rC_Z_g2, rC_Z_g2);
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
        copy_libsnark(proof.g_A.g, AG);
        copy_libsnark(proof.g_A.h, AH);
        copy_libsnark(proof.g_B.g, BG);
        copy_libsnark(proof.g_B.h, BH);
        copy_libsnark(proof.g_C.g, CG);
        copy_libsnark(proof.g_C.h, CH);
        copy_libsnark(proof.g_H, H);
        copy_libsnark(proof.g_K, K);
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
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// compare original and redesigned proof output
//

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_ProofCompare : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp PPT;
#else
    typedef libsnark::default_ec_pp PPT;
#endif

public:
    AutoTest_PPZK_ProofCompare(const AutoTestR1CS<SYS, Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const auto keypair
            = libsnark::r1cs_ppzksnark_generator<PPT>(
                m_constraintSystem.systemA());

#ifdef USE_OLD_LIBSNARK
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.witnessA());
#else
        const auto proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                keypair.pk,
                m_constraintSystem.inputA(),
                m_constraintSystem.witnessA());
#endif

        // proof from original code
        G1 AG, AH, BH, CG, CH, H, K;
        G2 BG;
        copy_libsnark(proof.g_A.g, AG);
        copy_libsnark(proof.g_A.h, AH);
        copy_libsnark(proof.g_B.g, BG);
        copy_libsnark(proof.g_B.h, BH);
        copy_libsnark(proof.g_C.g, CG);
        copy_libsnark(proof.g_C.h, CH);
        copy_libsnark(proof.g_H, H);
        copy_libsnark(proof.g_K, K);
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
        copy_libsnark(keypair.pk.A_query, A_query);
        copy_libsnark(keypair.pk.B_query, B_query);
        copy_libsnark(keypair.pk.C_query, C_query);
        copy_libsnark(keypair.pk.H_query, H_query);
        copy_libsnark(keypair.pk.K_query, K_query);
        const PPZK_ProvingKey<PAIRING> pkB(A_query,
                                           B_query,
                                           C_query,
                                           H_query,
                                           K_query);

        // proof from redesigned code
        const PPZK_Proof<PAIRING> proofFromRedesign(
            m_constraintSystem.systemB(),
            m_constraintSystem.numCircuitInputs(),
            pkB,
            m_constraintSystem.witnessB(),
            PPZK_ProofRandomness<typename PAIRING::Fr>(0));

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
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// verification and proof use redesigned code
//

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_Proof : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp PPT;
#else
    typedef libsnark::default_ec_pp PPT;
#endif

public:
    AutoTest_PPZK_Proof(const AutoTestR1CS<SYS, Fr, U>& cs)
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
#ifdef USE_OLD_LIBSNARK
        copy_libsnark(keypair.pk.A_query, A_query);
        copy_libsnark(keypair.pk.B_query, B_query);
        copy_libsnark(keypair.pk.C_query, C_query);
        copy_libsnark(keypair.pk.H_query, H_query);
        copy_libsnark(keypair.pk.K_query, K_query);
#else
        // new and old libsnark have different query vector formats
        // new libsnark appends inhomogeneous Z values to the back
        // old libsnark prepends inhomogeneous Z values to the front
        // convert new format to old format which snarklib expects

        // A
        {
            const auto Z = keypair.pk.A_query.values.back();
            const auto len = keypair.pk.A_query.size();
            A_query.reserve(len + 2);
            G1 tmpG, tmpH;
            copy_libsnark(Z.g, tmpG);
            copy_libsnark(Z.h, tmpH);
            A_query.pushBack(0, Pairing<G1, G1>(tmpG, tmpH));
            A_query.pushBack(1, Pairing<G1, G1>::zero());
            A_query.pushBack(2, Pairing<G1, G1>::zero());
            for (std::size_t i = 0; i < len - 1; ++i) {
                copy_libsnark(keypair.pk.A_query.values[i].g, tmpG);
                copy_libsnark(keypair.pk.A_query.values[i].h, tmpH);

                A_query.pushBack(
                    keypair.pk.A_query.indices[i] + 3,
                    Pairing<G1, G1>(tmpG, tmpH));
            }
        }

        // B
        {
            const auto Z = keypair.pk.B_query.values.back();
            const auto len = keypair.pk.B_query.size();
            B_query.reserve(len + 2);
            G2 tmpG;
            G1 tmpH;
            copy_libsnark(Z.g, tmpG);
            copy_libsnark(Z.h, tmpH);
            B_query.pushBack(0, Pairing<G2, G1>::zero());
            B_query.pushBack(1, Pairing<G2, G1>(tmpG, tmpH));
            B_query.pushBack(2, Pairing<G2, G1>::zero());
            for (std::size_t i = 0; i < len - 1; ++i) {
                copy_libsnark(keypair.pk.B_query.values[i].g, tmpG);
                copy_libsnark(keypair.pk.B_query.values[i].h, tmpH);

                B_query.pushBack(
                    keypair.pk.B_query.indices[i] + 3,
                    Pairing<G2, G1>(tmpG, tmpH));
            }
        }

        // C
        {
            const auto Z = keypair.pk.C_query.values.back();
            const auto len = keypair.pk.C_query.size();
            C_query.reserve(len + 2);
            G1 tmpG, tmpH;
            copy_libsnark(Z.g, tmpG);
            copy_libsnark(Z.h, tmpH);
            C_query.pushBack(0, Pairing<G1, G1>::zero());
            C_query.pushBack(1, Pairing<G1, G1>::zero());
            C_query.pushBack(2, Pairing<G1, G1>(tmpG, tmpH));
            for (std::size_t i = 0; i < len - 1; ++i) {
                copy_libsnark(keypair.pk.C_query.values[i].g, tmpG);
                copy_libsnark(keypair.pk.C_query.values[i].h, tmpH);

                C_query.pushBack(
                    keypair.pk.C_query.indices[i] + 3,
                    Pairing<G1, G1>(tmpG, tmpH));
            }
        }

        // H
        copy_libsnark(keypair.pk.H_query, H_query);

        // K
        {
            const auto lenK = keypair.pk.K_query.size();
            K_query.reserve(lenK);
            G1 tmpG;
            copy_libsnark(keypair.pk.K_query[lenK - 3], tmpG);
            K_query.emplace_back(tmpG);
            copy_libsnark(keypair.pk.K_query[lenK - 2], tmpG);
            K_query.emplace_back(tmpG);
            copy_libsnark(keypair.pk.K_query[lenK - 1], tmpG);
            K_query.emplace_back(tmpG);
            for (std::size_t i = 0; i < lenK - 3; ++i) {
                copy_libsnark(keypair.pk.K_query[i], tmpG);
                K_query.emplace_back(tmpG);
            }
        }
#endif
        const PPZK_ProvingKey<PAIRING> pkB(A_query,
                                           B_query,
                                           C_query,
                                           H_query,
                                           K_query);

        // encoded IC query
        G1 base;
#ifdef USE_OLD_LIBSNARK
        const size_t encSize = keypair.vk.encoded_IC_query->encoded_terms.size();
#else
        const size_t encSize = keypair.vk.encoded_IC_query.rest.size();
#endif
        std::vector<G1> encoded_terms(encSize);
#ifdef USE_OLD_LIBSNARK
        copy_libsnark(keypair.vk.encoded_IC_query->base, base);
#else
        copy_libsnark(keypair.vk.encoded_IC_query.first, base);
#endif
        for (std::size_t i = 0; i < encSize; ++i) {
#ifdef USE_OLD_LIBSNARK
            copy_libsnark(keypair.vk.encoded_IC_query->encoded_terms[i], encoded_terms[i]);
#else
            copy_libsnark(keypair.vk.encoded_IC_query.rest[i], encoded_terms[i]);
#endif
        }
        const PPZK_QueryIC<PAIRING> icqB(base, encoded_terms);

        // verification key
        G1 alphaB_g1, gamma_beta_g1;
        G2 alphaA_g2, alphaC_g2, gamma_g2, gamma_beta_g2, rC_Z_g2;
        copy_libsnark(keypair.vk.alphaA_g2, alphaA_g2);
        copy_libsnark(keypair.vk.alphaB_g1, alphaB_g1);
        copy_libsnark(keypair.vk.alphaC_g2, alphaC_g2);
        copy_libsnark(keypair.vk.gamma_g2, gamma_g2);
        copy_libsnark(keypair.vk.gamma_beta_g1, gamma_beta_g1);
        copy_libsnark(keypair.vk.gamma_beta_g2, gamma_beta_g2);
        copy_libsnark(keypair.vk.rC_Z_g2, rC_Z_g2);
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
                                         m_constraintSystem.numCircuitInputs(),
                                         pkB,
                                         m_constraintSystem.witnessB(),
                                         PPZK_ProofRandomness<typename PAIRING::Fr>(0));

        const auto ans
            = strongVerify(
                vkB,
                m_constraintSystem.inputB(),
                proofB);

        checkPass(ans);
    }

private:
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// end-to-end using snarklib
//

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_full_redesign : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp PPT;
#else
    typedef libsnark::default_ec_pp PPT;
#endif

public:
    AutoTest_PPZK_full_redesign(const AutoTestR1CS<SYS, Fr, U>& cs,
                                const bool failure_is_success = false)
        : AutoTest(cs),
          m_constraintSystem(cs),
          m_failureIsSuccess(failure_is_success)
    {}

    void runTest() {
        const PPZK_Keypair<PAIRING> keypair(m_constraintSystem.systemB(),
                                            m_constraintSystem.numCircuitInputs(),
                                            PPZK_LagrangePoint<Fr>(0),
                                            PPZK_BlindGreeks<Fr, Fr>(0));

        const PPZK_Proof<PAIRING> proofB(m_constraintSystem.systemB(),
                                         m_constraintSystem.numCircuitInputs(),
                                         keypair.pk(),
                                         m_constraintSystem.witnessB(),
                                         PPZK_ProofRandomness<typename PAIRING::Fr>(0));

        const auto ans = strongVerify(keypair.vk(),
                                      m_constraintSystem.inputB(),
                                      proofB);

        checkPass(m_failureIsSuccess ? !ans : ans);
    }

private:
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
    const bool m_failureIsSuccess;
};

} // namespace snarklib

#endif
