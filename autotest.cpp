#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <gmp.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include "AutoTest.hpp"
#include "AutoTest_BigInt.hpp"
#include "AutoTest_EC_Pairing.hpp"
#include "AutoTest_Field.hpp"
#include "AutoTest_Group.hpp"
#include "AutoTest_LagrangeFFT.hpp"
#include "AutoTest_Marshalling.hpp"
#include "AutoTest_MultiExp.hpp"
#include "AutoTest_Pairing.hpp"
#include "AutoTest_PPZK.hpp"
#include "AutoTest_QAP.hpp"
#include "AutoTest_R1CS.hpp"
#include "AutoTest_WindowExp.hpp"
#include "common/profiling.hpp"
#include "common/types.hpp"
#include "EC.hpp"
#include "Field.hpp"
#include "FpModel.hpp"
#include "Rank1DSL.hpp"

using namespace snarklib;
using namespace std;

// libsnark is compiled with a hardcoded elliptic curve. The library
// can only use the one elliptic curve set during compilation.
#ifdef CURVE_ALT_BN128
typedef BN128 CURVE;
#endif
#ifdef CURVE_EDWARDS
typedef Edwards CURVE;
#endif

// elliptic curve types
const mp_size_t NRQ = CURVE::q_limbs;
extern const auto MODULUS_R = CURVE::modulus_r(); // note: external linkage
extern const auto MODULUS_Q = CURVE::modulus_q(); // note: external linkage

typedef typename CURVE::Groups<NRQ, MODULUS_R, MODULUS_Q>::Fr Fr;
typedef typename CURVE::Groups<NRQ, MODULUS_R, MODULUS_Q>::Fq Fq;
typedef typename CURVE::Groups<NRQ, MODULUS_R, MODULUS_Q>::G1 G1;
typedef typename CURVE::Groups<NRQ, MODULUS_R, MODULUS_Q>::G2 G2;

typedef typename CURVE::Pairing<NRQ, MODULUS_R, MODULUS_Q> PAIRING;
typedef typename PAIRING::Fqe Fqe;
typedef typename PAIRING::Fqk Fqk;
typedef typename PAIRING::G1_precomp G1_precomp;
typedef typename PAIRING::G2_precomp G2_precomp;
typedef typename PAIRING::GT GT;

// initialize elliptic curve parameters
void initEC() {
    // the R and Q modulus should be about the same size for GMP
    assert(CURVE::r_limbs == CURVE::q_limbs);

    // snarklib:
    // critically important to initialize finite field and group parameters
    CURVE::Fields<NRQ, MODULUS_R>::initParams();
    CURVE::Fields<NRQ, MODULUS_Q>::initParams();
    CURVE::Groups<NRQ, MODULUS_R, MODULUS_Q>::initParams();

    // libsnark:
    // critically important to initialize finite fiels and group parameters
    libsnark::init_public_params<libsnark::default_pp>();

    // disable/reduce instrumentation messages
    libsnark::inhibit_profiling_counters = true;
}

random_device rd;

template <mp_size_t N>
void add_BigInt(AutoTestBattery& ATB)
{
    ATB.addTest(new AutoTest_BigIntDefaultConstructorZero<N>);
    ATB.addTest(new AutoTest_BigIntUnsignedLongConstructor<N>(0, 1, 100));

    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_BigIntUnsignedLongConstructor<N>(rd(), rd(), 100));
        ATB.addTest(new AutoTest_BigIntStringConstructor<N>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_BigIntEquality<N>(randomBase10(rd, N), randomBase10(rd, N)));
        ATB.addTest(new AutoTest_BigIntStreamOutput<N>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_BigIntNumBits<N>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_BigIntTestBits<N>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_BigIntFindwNAF<N>(rd() % 16, randomBase10(rd, N)));
    }
}

template <mp_size_t N, typename T, typename U>
void add_Field(AutoTestBattery& ATB)
{
    ATB.addTest(new AutoTest_FieldZeroStaysZero<T>);
    ATB.addTest(new AutoTest_FieldOneStaysOne<T>);

    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_FieldZeroAndOneIdentities<T>);
        ATB.addTest(new AutoTest_FieldAdd<T, U>);
        ATB.addTest(new AutoTest_FieldSub<T, U>);
        ATB.addTest(new AutoTest_FieldMul<T, U>);
        ATB.addTest(new AutoTest_FieldExp<N, T, U>);
        ATB.addTest(new AutoTest_FieldSquared<T, U>);
        ATB.addTest(new AutoTest_FieldInverse<T, U>);
    }
}

template <typename T, typename U>
void add_Field_sqrt(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        // defined for only: Fp, Fp2, Fp3
        ATB.addTest(new AutoTest_FieldSqrt<T, U>);
    }
}

template <typename T, typename U>
void add_Field_Frobenius_map(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        // defined for everything except Fp
        ATB.addTest(new AutoTest_FieldFrobeniusMap<T, U>(rd()));
    }
}

template <mp_size_t N, typename T, typename U>
void add_Field_cyclotomic_exp(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        // defined for Fp32 and Fp232, the Fqk/GT fields for BN128 and Edwards
        ATB.addTest(new AutoTest_FieldCyclotomicExp<N, T, U>);
    }
}

#ifdef CURVE_ALT_BN128
void add_Field_mul_by_024(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(
            new AutoTest_FieldMulBy024<Fqk,
                                       Fqe,
                                       libsnark::Fqk<libsnark::default_pp>,
                                       libsnark::Fqe<libsnark::default_pp>>);
    }
}
#endif

template <mp_size_t N, typename T, typename U>
void add_Group(AutoTestBattery& ATB)
{
    ATB.addTest(new AutoTest_GroupZeroStaysZero<T>);

    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_GroupZeroIdentity<T>);
        ATB.addTest(new AutoTest_GroupOneIdentity<T>);
        ATB.addTest(new AutoTest_GroupAdd<N, T, U>(randomBase10(rd, N), randomBase10(rd, N)));
        ATB.addTest(new AutoTest_GroupSub<N, T, U>(randomBase10(rd, N), randomBase10(rd, N)));
        ATB.addTest(new AutoTest_GroupMul<N, T, U>(randomBase10(rd, N), randomBase10(rd, N)));
        ATB.addTest(new AutoTest_GroupDbl<N, T, U>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_GroupSpecialWellFormed<N, T, U>(randomBase10(rd, N)));
    }
}

template <mp_size_t N, typename PAIRING, typename UG1, typename UG2, typename UGT>
void add_EC_Pairing(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_EC_PairingPrecompG1<N, PAIRING, UG1>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_EC_PairingPrecompG2<N, PAIRING, UG2>(randomBase10(rd, N)));
        ATB.addTest(new AutoTest_EC_PairingDoublingStepForFlippedMillerLoop<N, PAIRING, UG2>(
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_EC_PairingAteMillerLoop<N, PAIRING, UG1, UG2>(
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_EC_PairingAteDoubleMillerLoop<N, PAIRING, UG1, UG2>(
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_EC_PairingFinalExponentiation<N, PAIRING, UGT>);
    }
}

template <mp_size_t N, typename T, typename F, typename U, typename G>
void add_MultiExp(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_MultiExp_wnafExp<N, T, U>(
                        uniformBase10(0, 1000000),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_MultiExp_multiExp<N, T, F, U, G>(rd() % 100));
        ATB.addTest(new AutoTest_MultiExp_multiExp01<N, T, F, U, G>(rd() % 100));
    }
}

template <typename T, typename U>
void add_LagrangeFFT(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_LagrangeFFT_FFT<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_iFFT<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_cosetFFT<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_icosetFFT<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_lagrange_coeffs<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_get_element<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_compute_Z<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_add_poly_Z<T, U>(2 + rd() % 100));
        ATB.addTest(new AutoTest_LagrangeFFT_divide_by_Z_on_coset<T, U>(2 + rd() % 100));
    }
}

template <typename T, typename F, typename U, typename G>
void add_WindowExp(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_WindowExpSize<T, U>(1 + rd() % 100));
        ATB.addTest(new AutoTest_WindowExp_exp<T, F, U, G>(1 + rd() % 100));
        ATB.addTest(new AutoTest_WindowExp_batchExp<T, F, U, G>(
                        1 + rd() % 100,
                        1 + rd() % 10));
    }

    for (size_t i = 0; i < 2; ++i) {
        ATB.addTest(new AutoTest_WindowExp_expMapReduce<T, F>(1 + rd() % 100));
        ATB.addTest(new AutoTest_WindowExp_batchExpMapReduce1<T, F>(
                        1 + rd() % 100,
                        1 + rd() % 10));
        ATB.addTest(new AutoTest_WindowExp_batchExpMapReduce2<T, F>(
                        1 + rd() % 100,
                        1 + rd() % 10));
    }
}

template <mp_size_t N,
          typename TG, typename TH, typename TF,
          typename UG, typename UH, typename UF>
void add_Pairing(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 10; ++i) {
        ATB.addTest(new AutoTest_PairingAdd<N, TG, TH, UG, UH>(
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_PairingMul<N, TG, TH, TF, UG, UH, UF>(
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_PairingFastAddSpecial<N, TG, TH, UG, UH>(
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_PairingBatchSpecial<N, TG, TH, UG, UH>(1 + rd() % 100));
        ATB.addTest(new AutoTest_Pairing_wnafExp<N, TG, TH, UG, UH>(
                        uniformBase10(0, 1000000),
                        randomBase10(rd, N),
                        randomBase10(rd, N)));
        ATB.addTest(new AutoTest_Pairing_batchExp<N, TG, TH, TF, UG, UH, UF>(
                        1 + rd() % 100,
                        1 + rd() % 10));
        ATB.addTest(new AutoTest_Pairing_multiExp01<N, TG, TH, TF, UG, UH, UF>(
                        1 + rd() % 100));
    }

    for (size_t i = 0; i < 2; ++i) {
        ATB.addTest(new AutoTest_Pairing_batchExpMapReduce1<TG, TH, TF>(
                        1 + rd() % 100,
                        1 + rd() % 10));
        ATB.addTest(new AutoTest_Pairing_batchExpMapReduce2<TG, TH, TF>(
                        1 + rd() % 100,
                        1 + rd() % 10));
    }
}

template <typename T, typename U>
void add_QAP(AutoTestBattery& ATB)
{
    for (const auto x_IC : { false, true }) {
        for (const auto y_IC : { false, true}) {
            const std::vector<AutoTestR1CS<T, U>> csvec = {
                AutoTestR1CS_AND<T, U>(x_IC, y_IC),
                AutoTestR1CS_OR<T, U>(x_IC, y_IC),
                AutoTestR1CS_XOR<T, U>(x_IC, y_IC),
                AutoTestR1CS_CMPLMNT<T, U>(x_IC)
            };

            for (size_t i = 0; i < 2; ++i) {
                for (const auto& cs : csvec) {
                    ATB.addTest(new AutoTest_QAP_ABCH_instance_map<T, U>(cs));
                    ATB.addTest(new AutoTest_QAP_Witness_map<T, U>(cs));
                }
            }
        }
    }
}

template <typename PAIRING, typename T, typename U>
void add_PPZK(AutoTestBattery& ATB)
{
    for (const auto x_IC : { false, true }) {
        for (const auto y_IC : { false, true}) {
            const std::vector<AutoTestR1CS<T, U>> csvec = {
                AutoTestR1CS_AND<T, U>(x_IC, y_IC),
                AutoTestR1CS_OR<T, U>(x_IC, y_IC),
                AutoTestR1CS_XOR<T, U>(x_IC, y_IC),
                AutoTestR1CS_CMPLMNT<T, U>(x_IC)
            };

            for (size_t i = 0; i < 2; ++i) {
                for (const auto& cs : csvec) {
                    ATB.addTest(new AutoTest_PPZK_libsnark_only<T, U>(cs));
                    ATB.addTest(new AutoTest_PPZK_strongVerify<PAIRING, U>(cs));
                    ATB.addTest(new AutoTest_PPZK_ProofCompare<PAIRING, U>(cs));
                    ATB.addTest(new AutoTest_PPZK_Proof<PAIRING, U>(cs));
                    ATB.addTest(new AutoTest_PPZK_full_redesign<PAIRING, U>(cs));
                }
            }
        }
    }
}

template <typename GA, typename GB, mp_size_t N, typename F, typename PAIRING>
void add_Marshalling(AutoTestBattery& ATB)
{
    for (size_t i = 0; i < 2; ++i) {
        ATB.addTest(new AutoTest_Marshal_SparseVectorPairing<GA, GB>(rd() % 100, rd() % 10));
        ATB.addTest(new AutoTest_Marshal_BFG<BigInt<N>>);
        ATB.addTest(new AutoTest_Marshal_BFG<F>);
        ATB.addTest(new AutoTest_Marshal_BFG<GA>);
        ATB.addTest(new AutoTest_Marshal_BFG<GB>);
        ATB.addTest(new AutoTest_Marshal_Pairing<GA, GB>);
        ATB.addTest(new AutoTest_Marshal_ProvingKey<PAIRING>(rd() % 100, rd() % 10));
        ATB.addTest(new AutoTest_Marshal_QueryIC<PAIRING>(rd() % 100));
        ATB.addTest(new AutoTest_Marshal_VerificationKey<PAIRING>(rd() % 100));
        ATB.addTest(new AutoTest_Marshal_Keypair<PAIRING>(rd() % 100, rd() % 10));
        ATB.addTest(new AutoTest_Marshal_Proof<PAIRING>);
        ATB.addTest(new AutoTest_Marshal_R1Witness<F>(rd() % 100));
    }
}

void printUsage(const char* exeName) {
    cout << "run all tests:  " << exeName << " -a" << endl
         << "specified test: " << exeName << " -i testnumber" << endl;

    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    // command line switches
    bool ok = false;
    size_t testNumber = -1;
    int opt;
    while (-1 != (opt = getopt(argc, argv, "ai:"))) {
        switch (opt) {
        case ('a') :
            ok = true;
            break;
        case ('i') :
            {
                stringstream ss(optarg);
                ss >> testNumber;
                ok = !!ss;
            }
            break;
        }
    }

    if (!ok) printUsage(argv[0]);

    // initialize snarklib and libsnark
    initEC();

    typedef typename libsnark::Fr<libsnark::default_pp> libsnark_Fr;

    typedef typename libsnark::Fq<libsnark::default_pp> libsnark_Fq;
    typedef typename libsnark::Fqe<libsnark::default_pp> libsnark_Fqe;
    typedef typename libsnark::Fqk<libsnark::default_pp> libsnark_Fqk;

    typedef typename libsnark::G1<libsnark::default_pp> libsnark_G1;
    typedef typename libsnark::G2<libsnark::default_pp> libsnark_G2;

    AutoTestBattery ATB;

    // big integers
    add_BigInt<1>(ATB);
    add_BigInt<2>(ATB);
    add_BigInt<3>(ATB);
    add_BigInt<4>(ATB);
    add_BigInt<5>(ATB);
    add_BigInt<6>(ATB);

    // algebraic fields
    add_Field<NRQ, Fr, libsnark_Fr>(ATB);
    add_Field<NRQ, Fq, libsnark_Fq>(ATB);
    add_Field<NRQ, Fqe, libsnark_Fqe>(ATB);
    add_Field<NRQ, Fqk, libsnark_Fqk>(ATB);
    add_Field_sqrt<Fr, libsnark_Fr>(ATB);
    add_Field_sqrt<Fq, libsnark_Fq>(ATB);
    add_Field_sqrt<Fqe, libsnark_Fqe>(ATB);
    add_Field_Frobenius_map<Fqe, libsnark_Fqe>(ATB);
    add_Field_Frobenius_map<Fqk, libsnark_Fqk>(ATB);
    add_Field_cyclotomic_exp<NRQ, Fqk, libsnark_Fqk>(ATB);
#ifdef CURVE_ALT_BN128
    add_Field_mul_by_024(ATB);
#endif

    // algebraic groups
    add_Group<NRQ, G1, libsnark_G1>(ATB);
    add_Group<NRQ, G2, libsnark_G2>(ATB);

    // pairings
    add_EC_Pairing<NRQ, PAIRING, libsnark_G1, libsnark_G2, libsnark_Fqk>(ATB);

    // multiple exponentiation
    add_MultiExp<NRQ, G1, Fr, libsnark_G1, libsnark_Fr>(ATB);
    add_MultiExp<NRQ, G2, Fr, libsnark_G2, libsnark_Fr>(ATB);

    // Lagrange FFT
    add_LagrangeFFT<Fr, libsnark_Fr>(ATB);

    // window table
    add_WindowExp<G1, Fr, libsnark_G1, libsnark_Fr>(ATB);
    add_WindowExp<G2, Fr, libsnark_G2, libsnark_Fr>(ATB);

    // paired groups
    add_Pairing<NRQ, G1, G1, Fr, libsnark_G1, libsnark_G1, libsnark_Fr>(ATB);
    add_Pairing<NRQ, G2, G1, Fr, libsnark_G2, libsnark_G1, libsnark_Fr>(ATB);

    // quadratic arithmetic program
    add_QAP<Fr, libsnark_Fr>(ATB);

    // pre-processed zero knowledge proof
    add_PPZK<PAIRING, Fr, libsnark_Fr>(ATB);

    // marshalling
    add_Marshalling<G1, G1, 1, Fr, PAIRING>(ATB);
    add_Marshalling<G1, G2, 2, Fr, PAIRING>(ATB);
    add_Marshalling<G2, G1, 3, Fq, PAIRING>(ATB);
    add_Marshalling<G2, G2, 4, Fq, PAIRING>(ATB);

    if (-1 == testNumber) {
        // run all tests
        if (ATB.runTest())
            cout << ATB.testCount() << " tests passed" << endl;
        else
            ATB.testLog(cout);

    } else {
        // run specified test only
        ATB.runTest(testNumber);
        ATB.testLog(cout, testNumber);
    }

    exit(EXIT_SUCCESS);
}
