#ifndef _SNARKLIB_AUTOTEST_PPZK_HPP_
#define _SNARKLIB_AUTOTEST_PPZK_HPP_

#include <cstdint>
#include <vector>

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
    typedef LIBSNARK_PPT PPT;

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

    typedef LIBSNARK_PPT PPT;

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

        // verification key
        PPZK_VerificationKey<PAIRING> vkB;
        copy_libsnark(keypair.vk, vkB);

        // proof
        PPZK_Proof<PAIRING> proofB;
        copy_libsnark(proof, proofB);

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

    typedef LIBSNARK_PPT PPT;

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
        PPZK_Proof<PAIRING> proofFromOriginal;
        copy_libsnark(proof, proofFromOriginal);

        // proving key
        PPZK_ProvingKey<PAIRING> pkB;
        copy_libsnark(keypair.pk, pkB);

        // proof from redesigned code
        const PPZK_Proof<PAIRING> proofFromRedesign(
            m_constraintSystem.systemB(),
            m_constraintSystem.numCircuitInputs(),
            pkB,
            m_constraintSystem.witnessB(),
            PPZK_ProofRandomness<Fr>(0));

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

    typedef LIBSNARK_PPT PPT;

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
        PPZK_ProvingKey<PAIRING> pkB;
        copy_libsnark(keypair.pk, pkB);

        // verification key
        PPZK_VerificationKey<PAIRING> vkB;
        copy_libsnark(keypair.vk, vkB);

        // proof
        const PPZK_Proof<PAIRING> proofB(m_constraintSystem.systemB(),
                                         m_constraintSystem.numCircuitInputs(),
                                         pkB,
                                         m_constraintSystem.witnessB(),
                                         PPZK_ProofRandomness<Fr>(0));

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
                                         PPZK_ProofRandomness<Fr>(0));

        const auto ans = strongVerify(keypair.vk(),
                                      m_constraintSystem.inputB(),
                                      proofB);

        checkPass(m_failureIsSuccess ? !ans : ans);
    }

private:
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
    const bool m_failureIsSuccess;
};

////////////////////////////////////////////////////////////////////////////////
// verification uses libsnark code
//

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_strongVerify_libsnark : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef LIBSNARK_PPT PPT;
    typedef LIBSNARK_FR FR;

public:
    AutoTest_PPZK_strongVerify_libsnark(const AutoTestR1CS<SYS, Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
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
                                         PPZK_ProofRandomness<Fr>(0));

        libsnark::r1cs_ppzksnark_verification_key<PPT> libsnark_vk;
        copy_libsnark(keypair.vk(), libsnark_vk);

        std::vector<FR> libsnark_input;
        copy_libsnark(m_constraintSystem.inputB(), libsnark_input);

        libsnark::r1cs_ppzksnark_proof<PPT> libsnark_proof;
        copy_libsnark(proofB, libsnark_proof);

        const auto ans_online
            = libsnark::r1cs_ppzksnark_verifier_strong_IC<PPT>(
                libsnark_vk,
                libsnark_input,
                libsnark_proof);

        checkPass(ans_online);
    }

private:
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
};

////////////////////////////////////////////////////////////////////////////////
// verification and proof use libsnark code
//

template <template <typename> class SYS, typename PAIRING, typename U>
class AutoTest_PPZK_Proof_libsnark : public AutoTest
{
    typedef typename PAIRING::Fr Fr;
    typedef LIBSNARK_PPT PPT;
    typedef LIBSNARK_FR FR;

public:
    AutoTest_PPZK_Proof_libsnark(const AutoTestR1CS<SYS, Fr, U>& cs)
        : AutoTest(cs),
          m_constraintSystem(cs)
    {}

    void runTest() {
        const PPZK_Keypair<PAIRING> keypair(m_constraintSystem.systemB(),
                                            m_constraintSystem.numCircuitInputs(),
                                            PPZK_LagrangePoint<Fr>(0),
                                            PPZK_BlindGreeks<Fr, Fr>(0));

        libsnark::r1cs_ppzksnark_proving_key<PPT> libsnark_pk;
        copy_libsnark(keypair.pk(), m_constraintSystem.systemA(), libsnark_pk);

        libsnark::r1cs_ppzksnark_verification_key<PPT> libsnark_vk;
        copy_libsnark(keypair.vk(), libsnark_vk);

        std::vector<FR> libsnark_input;
        copy_libsnark(m_constraintSystem.inputB(), libsnark_input);

        std::vector<FR> libsnark_witness;
#ifdef USE_OLD_LIBSNARK
        copy_libsnark(m_constraintSystem.witnessB(), libsnark_witness);
#else
        const std::size_t numInputs = libsnark_input.size();
        copy_libsnark(m_constraintSystem.witnessB(), libsnark_witness, numInputs);
#endif

#ifdef USE_OLD_LIBSNARK
        const auto libsnark_proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                libsnark_pk,
                libsnark_witness);
#else
        const auto libsnark_proof
            = libsnark::r1cs_ppzksnark_prover<PPT>(
                libsnark_pk,
                libsnark_input,
                libsnark_witness);
#endif

        const auto ans_online
            = libsnark::r1cs_ppzksnark_verifier_strong_IC<PPT>(
                libsnark_vk,
                libsnark_input,
                libsnark_proof);

        checkPass(ans_online);
    }

private:
    const AutoTestR1CS<SYS, Fr, U> m_constraintSystem;
};

} // namespace snarklib

#endif
