#ifndef _SNARKLIB_FOREIGN_LIB_HPP_
#define _SNARKLIB_FOREIGN_LIB_HPP_

#include <algorithm>
#include <cstdint>
#include <gmp.h>
#include <sstream>
#include <vector>

#include /*libsnark*/ "algebra/curves/alt_bn128/alt_bn128_g1.hpp"
#include /*libsnark*/ "algebra/curves/alt_bn128/alt_bn128_g2.hpp"
#include /*libsnark*/ "algebra/curves/edwards/edwards_g1.hpp"
#include /*libsnark*/ "algebra/curves/edwards/edwards_g2.hpp"
#include /*libsnark*/ "algebra/fields/bigint.hpp"
#include /*libsnark*/ "algebra/fields/fp.hpp"
#include /*libsnark*/ "algebra/fields/fp2.hpp"
#include /*libsnark*/ "algebra/fields/fp3.hpp"
#include /*libsnark*/ "algebra/fields/fp6_2over3.hpp"
#include /*libsnark*/ "algebra/fields/fp6_3over2.hpp"
#include /*libsnark*/ "algebra/fields/fp12_2over3over2.hpp"

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "common/types.hpp"
#include /*libsnark*/ "r1cs_ppzksnark/r1cs_ppzksnark/hpp"
#include /*libsnark*/ "encoding/knowledge_commitment.hpp"
#include /*libsnark*/ "r1cs/r1cs.hpp"
#else
#include /*libsnark*/ "common/default_types/ec_pp.hpp"
#include /*libsnark*/ "zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp"
#include /*libsnark*/ "algebra/knowledge_commitment/knowledge_commitment.hpp"
#include /*libsnark*/ "relations/constraint_satisfaction_problems/r1cs/r1cs.hpp"
#endif

#include "snarklib/AuxSTL.hpp"
#include "snarklib/BigInt.hpp"
#include "snarklib/EC.hpp"
#include "snarklib/Field.hpp"
#include "snarklib/FpModel.hpp"
#include "snarklib/Group.hpp"
#include "snarklib/Pairing.hpp"
#include "snarklib/Rank1DSL.hpp"
#include "snarklib/PPZK_keypair.hpp"
#include "snarklib/PPZK_query.hpp"
#include "snarklib/PPZK_proof.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// elliptic curve typedefs
//

#ifdef USE_OLD_LIBSNARK
    typedef libsnark::default_pp LIBSNARK_PPT;
#else
    typedef libsnark::default_ec_pp LIBSNARK_PPT;
#endif

typedef libsnark::Fr<LIBSNARK_PPT> LIBSNARK_FR;
typedef libsnark::G1<LIBSNARK_PPT> LIBSNARK_G1;
typedef libsnark::G2<LIBSNARK_PPT> LIBSNARK_G2;

////////////////////////////////////////////////////////////////////////////////
// equality of libsnark and snarklib data structures
//

//
// bigint == BigInt
//

template <mp_size_t N>
bool equal_libsnark(
    const libsnark::bigint<N>& a,
    const BigInt<N>& b)
{
    for (std::size_t i = 0; i < N; ++i) {
        if (a.data[i] != b.data()[i])
            return false;
    }

    return true;
}

template <mp_size_t N>
bool equal_libsnark(
    const BigInt<N>& b,
    const libsnark::bigint<N>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp_model == FpModel
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp_model<N, MODULUS_A>& a,
    const FpModel<N, MODULUS_B>& b)
{
    return equal_libsnark(a.as_bigint(), b.asBigInt());
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const FpModel<N, MODULUS_B>& b,
    const libsnark::Fp_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp_model == Field<FpModel>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp_model<N, MODULUS_A>& a,
    const Field<FpModel<N, MODULUS_B>>& b)
{
    return equal_libsnark(a, b[0]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<FpModel<N, MODULUS_B>>& b,
    const libsnark::Fp_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp2_model == Field<FpModel, 2>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp2_model<N, MODULUS_A>& a,
    const Field<FpModel<N, MODULUS_B>, 2>& b)
{
    return
        equal_libsnark(a.c0, b[0]) &&
        equal_libsnark(a.c1, b[1]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<FpModel<N, MODULUS_B>, 2>& b,
    const libsnark::Fp2_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp3_model == Field<FpModel, 3>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp3_model<N, MODULUS_A>& a,
    const Field<FpModel<N, MODULUS_B>, 3>& b)
{
    return
        equal_libsnark(a.c0, b[0]) &&
        equal_libsnark(a.c1, b[1]) &&
        equal_libsnark(a.c2, b[2]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<FpModel<N, MODULUS_B>, 3>& b,
    const libsnark::Fp3_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp6_2over3_model == Field<Field<FpModel, 3>, 2>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp6_2over3_model<N, MODULUS_A>& a,
    const Field<Field<FpModel<N, MODULUS_B>, 3>, 2>& b)
{
    return
        equal_libsnark(a.c0, b[0]) &&
        equal_libsnark(a.c1, b[1]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<Field<FpModel<N, MODULUS_B>, 3>, 2>& b,
    const libsnark::Fp6_2over3_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp6_3over2_model == Field<Field<FpModel, 2>, 3>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp6_3over2_model<N, MODULUS_A>& a,
    const Field<Field<FpModel<N, MODULUS_B>, 2>, 3>& b)
{
    return
        equal_libsnark(a.c0, b[0]) &&
        equal_libsnark(a.c1, b[1]) &&
        equal_libsnark(a.c2, b[2]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<Field<FpModel<N, MODULUS_B>, 2>, 3>& b,
    const libsnark::Fp6_3over2_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// Fp12_2over3over2_model == Field<Field<Field<FpModel, 2>, 3>, 2>
//

template <mp_size_t N,
          const libsnark::bigint<N>& MODULUS_A,
          const BigInt<N>& MODULUS_B>
bool equal_libsnark(
    const libsnark::Fp12_2over3over2_model<N, MODULUS_A>& a,
    const Field<Field<Field<FpModel<N, MODULUS_B>, 2>, 3>, 2>& b)
{
    return
        equal_libsnark(a.c0, b[0]) &&
        equal_libsnark(a.c1, b[1]);
}

template <mp_size_t N,
          const BigInt<N>& MODULUS_B,
          const libsnark::bigint<N>& MODULUS_A>
bool equal_libsnark(
    const Field<Field<Field<FpModel<N, MODULUS_B>, 2>, 3>, 2>& b,
    const libsnark::Fp12_2over3over2_model<N, MODULUS_A>& a)
{
    return equal_libsnark(a, b);
}

//
// alt_bn128_G1 == BN128::Groups<>::G1
//

template <typename GROUP>
bool equal_libsnark(
    const libsnark::alt_bn128_G1& a,
    const GROUP& b)
{
    return
        equal_libsnark(a.X, b.x()) &&
        equal_libsnark(a.Y, b.y()) &&
        equal_libsnark(a.Z, b.z());
}

template <typename GROUP>
bool equal_libsnark(
    const GROUP& b,
    const libsnark::alt_bn128_G1& a)
{
    return equal_libsnark(a, b);
}

//
// alt_bn128_G2 == BN128::Groups<>::G2
//

template <typename GROUP>
bool equal_libsnark(
    const libsnark::alt_bn128_G2& a,
    const GROUP& b)
{
    return
        equal_libsnark(a.X, b.x()) &&
        equal_libsnark(a.Y, b.y()) &&
        equal_libsnark(a.Z, b.z());
}

template <typename GROUP>
bool equal_libsnark(
    const GROUP& b,
    const libsnark::alt_bn128_G2& a)
{
    return equal_libsnark(a, b);
}

//
// edwards_G1 == Edwards::Groups<>::G1
//

template <typename GROUP>
bool equal_libsnark(
    const libsnark::edwards_G1& a,
    const GROUP& b)
{
    return
        equal_libsnark(a.X, b.x()) &&
        equal_libsnark(a.Y, b.y()) &&
        equal_libsnark(a.Z, b.z());
}

template <typename GROUP>
bool equal_libsnark(
    const GROUP& b,
    const libsnark::edwards_G1& a)
{
    return equal_libsnark(a, b);
}

//
// edwards_G2 == Edwards::Groups<>::G2
//

template <typename GROUP>
bool equal_libsnark(
    const libsnark::edwards_G2& a,
    const GROUP& b)
{
    return
        equal_libsnark(a.X, b.x()) &&
        equal_libsnark(a.Y, b.y()) &&
        equal_libsnark(a.Z, b.z());
}

template <typename GROUP>
bool equal_libsnark(
    const GROUP& b,
    const libsnark::edwards_G2& a)
{
    return equal_libsnark(a, b);
}

//
// knowledge_commitment<> == Pairing<>
//

template <typename UG,
          typename UH,
          typename TG,
          typename TH>
bool equal_libsnark(
    const libsnark::knowledge_commitment<UG, UH>& a,
    const Pairing<TG, TH>& b)
{
    return
        equal_libsnark(a.g, b.G()) &&
        equal_libsnark(a.h, b.H());
}

template <typename TG,
          typename TH,
          typename UG,
          typename UH>
bool equal_libsnark(
    const Pairing<TG, TH>& b,
    const libsnark::knowledge_commitment<UG, UH>& a)
{
    return equal_libsnark(a, b);
}

//
// vector<> == vector<>
//

template <typename T,
          typename U>
bool equal_libsnark(
    const std::vector<T>& a,
    const std::vector<U>& b)
{
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (! equal_libsnark(a[i], b[i]))
            return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
// copy between libsnark and snarklib data structures
//

//
// bigint -> BigInt
//

template <mp_size_t N>
void copy_libsnark(
    const libsnark::bigint<N>& a,
    BigInt<N>& b)
{
    std::stringstream ssA;
    ssA << a;
    std::stringstream ssB(ssA.str());
    ssB >> b;
}

//
// BigInt -> bigint
//

template <mp_size_t N>
void copy_libsnark(
    const BigInt<N>& a,
    libsnark::bigint<N>& b)
{
    std::stringstream ssA;
    ssA << a;
    std::stringstream ssB(ssA.str());
    ssB >> b;
}

#define COMMA ,

//
// from libsnark FpX_model to snarklib Field<>
// 

#define DEFN_COPY_DATA_M2F(A, B)                        \
    template <mp_size_t N,                              \
              const libsnark::bigint<N>& MA,            \
              const BigInt<N>& MB>                      \
    void copy_libsnark(                                 \
        const libsnark:: A & a,                         \
        B & b)                                          \
    {                                                   \
        std::stringstream ssA;                          \
        ssA << a;                                       \
        std::stringstream ssB(ssA.str());               \
        ssB >> b;                                       \
    }

DEFN_COPY_DATA_M2F(Fp_model<N COMMA MA>, Field<FpModel<N COMMA MB>>)
DEFN_COPY_DATA_M2F(Fp2_model<N COMMA MA>, Field<FpModel<N COMMA MB> COMMA 2>)
DEFN_COPY_DATA_M2F(Fp3_model<N COMMA MA>, Field<FpModel<N COMMA MB> COMMA 3>)

// special case for libsnark Fp32 model to snarklib Field<>
#ifdef USE_OLD_LIBSNARK
template <mp_size_t N,
          const libsnark::bigint<N>& MA,
          const BigInt<N>& MB>
void copy_libsnark(
    const libsnark::Fp6_2over3_model<N, MA>& a,
    Field<Field<FpModel<N, MB>, 3>, 2>& b)
{
    std::stringstream ssA;
#ifdef CURVE_EDWARDS
    // Fp6_2over3_model has no overloaded operators for stream
    // insertion and extraction.
    ssA << a.c0 << " " << a.c1;
#else
    ssA << a;
#endif
    std::stringstream ssB(ssA.str());
    ssB >> b;
}
#else
DEFN_COPY_DATA_M2F(Fp6_2over3_model<N COMMA MA>, Field<Field<FpModel<N COMMA MB> COMMA 3> COMMA 2>)
#endif

DEFN_COPY_DATA_M2F(Fp6_3over2_model<N COMMA MA>, Field<Field<FpModel<N COMMA MB> COMMA 2> COMMA 3>)
DEFN_COPY_DATA_M2F(Fp12_2over3over2_model<N COMMA MA>, Field<Field<Field<FpModel<N COMMA MB> COMMA 2> COMMA 3> COMMA 2>)

#undef DEFN_COPY_DATA_M2F

//
// from snarklib Field<> to libsnark FpX_model
// 

#define DEFN_COPY_DATA_F2M(A, B)                        \
    template <mp_size_t N,                              \
              const BigInt<N>& MA,                      \
              const libsnark::bigint<N>& MB>            \
    void copy_libsnark(                                 \
        const A & a,                                    \
        libsnark:: B & b)                               \
    {                                                   \
        std::stringstream ssA;                          \
        ssA << a;                                       \
        std::stringstream ssB(ssA.str());               \
        ssB >> b;                                       \
    }

DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA>>, Fp_model<N COMMA MB>)
DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA> COMMA 2>, Fp2_model<N COMMA MB>)
DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA> COMMA 3>, Fp3_model<N COMMA MB>)

// special case for snarklib Field<> to libsnark Fp32 model
#ifdef USE_OLD_LIBSNARK
template <mp_size_t N,
          const BigInt<N>& MA,
          const libsnark::bigint<N>& MB>
void copy_libsnark(
    const Field<Field<FpModel<N, MA>, 3>, 2>& a,
    libsnark::Fp6_2over3_model<N, MB>& b)
{
    std::stringstream ssA;
    ssA << a;
    std::stringstream ssB(ssA.str());
#ifdef CURVE_EDWARDS
    // Fp6_2over3_model has no overloaded operators for stream
    // insertion and extraction.
    ssB >> b.c0 >> b.c1;
#else
    ssB >> b;
#endif
}
#else
DEFN_COPY_DATA_F2M(Field<Field<FpModel<N COMMA MA> COMMA 3> COMMA 2>, Fp6_2over3_model<N COMMA MB>)
#endif

DEFN_COPY_DATA_F2M(Field<Field<FpModel<N COMMA MA> COMMA 2> COMMA 3>, Fp6_3over2_model<N COMMA MB>)
DEFN_COPY_DATA_F2M(Field<Field<Field<FpModel<N COMMA MA> COMMA 2> COMMA 3> COMMA 2>, Fp12_2over3over2_model<N COMMA MB>)

#undef DEFN_COPY_DATA_F2M

//
// from libsnark elliptic curve group to snarklib Group<>
//

#define DEFN_COPY_DATA_M2G(A)                           \
    template <typename GROUP>                           \
    void copy_libsnark(                                 \
        const libsnark:: A & a,                         \
        GROUP& b)                                       \
    {                                                   \
        typename GROUP::BaseField tmp;                  \
        copy_libsnark(a.X, tmp);                        \
        b.x(tmp);                                       \
        copy_libsnark(a.Y, tmp);                        \
        b.y(tmp);                                       \
        copy_libsnark(a.Z, tmp);                        \
        b.z(tmp);                                       \
    }

DEFN_COPY_DATA_M2G(alt_bn128_G1)
DEFN_COPY_DATA_M2G(alt_bn128_G2)
DEFN_COPY_DATA_M2G(edwards_G1)
DEFN_COPY_DATA_M2G(edwards_G2)

#undef DEFN_COPY_DATA_M2G

//
// from snarklib Group<> to libsnark elliptic curve group
//

#define DEFN_COPY_DATA_G2M(B)                           \
    template <typename GROUP>                           \
    void copy_libsnark(                                 \
        const GROUP& a,                                 \
        libsnark:: B & b)                               \
    {                                                   \
        copy_libsnark(a.x(), b.X);                      \
        copy_libsnark(a.y(), b.Y);                      \
        copy_libsnark(a.z(), b.Z);                      \
    }

DEFN_COPY_DATA_G2M(alt_bn128_G1)
DEFN_COPY_DATA_G2M(alt_bn128_G2)
DEFN_COPY_DATA_G2M(edwards_G1)
DEFN_COPY_DATA_G2M(edwards_G2)

#undef DEFN_COPY_DATA_G2M

#undef COMMA

//
// from libsnark group vector to snarklib vector
//

template <typename G1>
void copy_libsnark(
    const libsnark::G1_vector<LIBSNARK_PPT>& a,
    std::vector<G1>& b,
    const std::size_t startIndex = 0,
    const std::size_t stopIndex = -1)
{
    const std::size_t L = std::min(a.size(), stopIndex);

    if (b.empty()) b.reserve(L - startIndex);

    G1 tmp;

    for (std::size_t i = startIndex; i < L; ++i) {
        copy_libsnark(a[i], tmp);
        b.emplace_back(tmp);
    }
}

//
// from snarklib vector to libsnark group vector
//

template <typename G1>
void copy_libsnark(
    const std::vector<G1>& a,
    libsnark::G1_vector<LIBSNARK_PPT>& b,
    const std::size_t startIndex = 0,
    const std::size_t stopIndex = -1)
{
    const std::size_t L = std::min(a.size(), stopIndex);

    if (b.empty()) b.reserve(L - startIndex);

    LIBSNARK_G1 tmp;

    for (std::size_t i = startIndex; i < L; ++i) {
        copy_libsnark(a[i], tmp);
        b.emplace_back(tmp);
    }
}

//
// from libsnark knowledge commitment to snarklib Pairing
//

template <typename UG,
          typename UH,
          typename TG,
          typename TH>
void copy_libsnark(
    const libsnark::knowledge_commitment<UG, UH>& a,
    Pairing<TG, TH>& b)
{
    TG tmpG;
    TH tmpH;

    copy_libsnark(a.g, tmpG);
    copy_libsnark(a.h, tmpH);

    b = Pairing<TG, TH>(tmpG, tmpH);
}

//
// from snarklib Pairing to libsnark knowledge commitment
//

template <typename TG,
          typename TH,
          typename UG,
          typename UH>
void copy_libsnark(
    const Pairing<TG, TH>& a,
    libsnark::knowledge_commitment<UG, UH>& b)
{
    UG tmpG;
    UH tmpH;

    copy_libsnark(a.G(), tmpG);
    copy_libsnark(a.H(), tmpH);

    b = libsnark::knowledge_commitment<UG, UH>(tmpG, tmpH);
}

//
// from libsnark paired groups vector to snarklib sparse vector
//

template <typename G1>
void copy_libsnark(
#ifdef USE_OLD_LIBSNARK
    const libsnark::G1G1_knowledge_commitment_vector<LIBSNARK_PPT>& a,
#else
    const libsnark::knowledge_commitment_vector<LIBSNARK_G1, LIBSNARK_G1>& a,
#endif
    SparseVector<Pairing<G1, G1>>& b,
    const std::size_t startIndex = 0,
    const std::size_t stopIndex = -1)
{
    const std::size_t
        L = std::min(a.values.size(), stopIndex),
        offset = b.size();

    if (b.empty()) b.reserve(L - startIndex);

    Pairing<G1, G1> tmp;

    for (std::size_t i = startIndex; i < L; ++i) {
        copy_libsnark(a.values[i], tmp);
        b.pushBack(a.indices[i] + offset, tmp);
    }
}

template <typename G2,
          typename G1>
void copy_libsnark(
#ifdef USE_OLD_LIBSNARK
    const libsnark::G2G1_knowledge_commitment_vector<LIBSNARK_PPT>& a,
#else
    const libsnark::knowledge_commitment_vector<LIBSNARK_G2, LIBSNARK_G1>& a,
#endif
    SparseVector<Pairing<G2, G1>>& b,
    const std::size_t startIndex = 0,
    const std::size_t stopIndex = -1)
{
    const std::size_t
        L = std::min(a.values.size(), stopIndex),
        offset = b.size();

    if (b.empty()) b.reserve(L - startIndex);

    Pairing<G2, G1> tmp;

    for (std::size_t i = startIndex; i < L; ++i) {
        copy_libsnark(a.values[i], tmp);
        b.pushBack(a.indices[i] + offset, tmp);
    }
}

//
// from snarklib sparse vector to libsnark paired groups vector
//

template <typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G1, G1>>& a,
#ifdef USE_OLD_LIBSNARK
    libsnark::G1G1_knowledge_commitment_vector<LIBSNARK_PPT>& b,
#else
    libsnark::knowledge_commitment_vector<LIBSNARK_G1, LIBSNARK_G1>& b,
#endif
    const std::size_t startElementIndex = 0,
    const std::size_t stopElementIndex = -1)
{
    if (b.indices.empty()) b.indices.reserve(a.size());
    if (b.values.empty()) b.values.reserve(a.size());

    std::size_t maxIndex = 0;
    libsnark::knowledge_commitment<LIBSNARK_G1, LIBSNARK_G1> tmp;

    for (std::size_t i = 0; i < a.size(); ++i) {
        const std::size_t idx = a.getIndex(i);

        if (idx >= startElementIndex && idx < stopElementIndex) {
            maxIndex = std::max(maxIndex, idx);

            const auto& p = a.getElement(i);
            copy_libsnark(p, tmp);

            b.indices.push_back(idx - startElementIndex);
            b.values.emplace_back(tmp);
        }
    }

#ifdef USE_OLD_LIBSNARK
    b.is_sparse = true;
    b.original_size = maxIndex + 1;
#else
    b.domain_size_ = maxIndex + 1;
#endif
}

template <typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G1, G1>>& a,
    const std::size_t index_Z,
    const Pairing<G1, G1>& value_Z,
#ifdef USE_OLD_LIBSNARK
    libsnark::G1G1_knowledge_commitment_vector<LIBSNARK_PPT>& b,
#else
    libsnark::knowledge_commitment_vector<LIBSNARK_G1, LIBSNARK_G1>& b,
#endif
    const std::size_t startElementIndex = 0,
    const std::size_t stopElementIndex = -1)
{
    if (b.indices.empty()) b.indices.reserve(a.size() + 1);
    if (b.values.empty()) b.values.reserve(a.size() + 1);

    std::size_t maxIndex = 0;
    libsnark::knowledge_commitment<LIBSNARK_G1, LIBSNARK_G1> tmp;

    for (std::size_t i = 0; i < a.size(); ++i) {
        const std::size_t idx = a.getIndex(i);

        if (idx >= startElementIndex && idx < stopElementIndex) {
            maxIndex = std::max(maxIndex, idx);

            const auto& p = a.getElement(i);
            copy_libsnark(p, tmp);

            b.indices.push_back(idx - startElementIndex);
            b.values.emplace_back(tmp);
        }
    }

    maxIndex = std::max(maxIndex, index_Z);
    b.indices.push_back(index_Z);
    copy_libsnark(value_Z, tmp);
    b.values.emplace_back(tmp);

#ifdef USE_OLD_LIBSNARK
    b.is_sparse = true;
    b.original_size = maxIndex + 1;
#else
    b.domain_size_ = maxIndex + 1;
#endif
}

template <typename G2,
          typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G2, G1>>& a,
#ifdef USE_OLD_LIBSNARK
    libsnark::G2G1_knowledge_commitment_vector<LIBSNARK_PPT>& b,
#else
    libsnark::knowledge_commitment_vector<LIBSNARK_G2, LIBSNARK_G1>& b,
#endif
    const std::size_t startElementIndex = 0,
    const std::size_t stopElementIndex = -1)
{
    if (b.indices.empty()) b.indices.reserve(a.size());
    if (b.values.empty()) b.values.reserve(a.size());

    std::size_t maxIndex = 0;
    libsnark::knowledge_commitment<LIBSNARK_G2, LIBSNARK_G1> tmp;

    for (std::size_t i = 0; i < a.size(); ++i) {
        const std::size_t idx = a.getIndex(i);

        if (idx >= startElementIndex && idx < stopElementIndex) {
            maxIndex = std::max(maxIndex, idx);

            const auto& p = a.getElement(i);
            copy_libsnark(p, tmp);

            b.indices.push_back(idx - startElementIndex);
            b.values.emplace_back(tmp);
        }
    }

#ifdef USE_OLD_LIBSNARK
    b.is_sparse = true;
    b.original_size = maxIndex + 1;
#else
    b.domain_size_ = maxIndex + 1;
#endif
}

template <typename G2,
          typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G2, G1>>& a,
    const std::size_t index_Z,
    const Pairing<G2, G1>& value_Z,
#ifdef USE_OLD_LIBSNARK
    libsnark::G2G1_knowledge_commitment_vector<LIBSNARK_PPT>& b,
#else
    libsnark::knowledge_commitment_vector<LIBSNARK_G2, LIBSNARK_G1>& b,
#endif
    const std::size_t startElementIndex = 0,
    const std::size_t stopElementIndex = -1)
{
    if (b.indices.empty()) b.indices.reserve(a.size() + 1);
    if (b.values.empty()) b.values.reserve(a.size() + 1);

    std::size_t maxIndex = 0;
    libsnark::knowledge_commitment<LIBSNARK_G2, LIBSNARK_G1> tmp;

    for (std::size_t i = 0; i < a.size(); ++i) {
        const std::size_t idx = a.getIndex(i);

        if (idx >= startElementIndex && idx < stopElementIndex) {
            maxIndex = std::max(maxIndex, idx);

            const auto& p = a.getElement(i);
            copy_libsnark(p, tmp);

            b.indices.push_back(idx - startElementIndex);
            b.values.emplace_back(tmp);
        }
    }

    maxIndex = std::max(maxIndex, index_Z);
    b.indices.push_back(index_Z);
    copy_libsnark(value_Z, tmp);
    b.values.emplace_back(tmp);

#ifdef USE_OLD_LIBSNARK
    b.is_sparse = true;
    b.original_size = maxIndex + 1;
#else
    b.domain_size_ = maxIndex + 1;
#endif
}

// 
// from snarklib to libsnark rank-1 linear combination
//

template <typename FR,
          typename LIBSNARK_FR>
void copy_libsnark(
    const R1Combination<FR>& a,
    libsnark::linear_combination<LIBSNARK_FR>& b)
{
    for (const auto& t : a.terms()) {
        LIBSNARK_FR tmp;
        copy_libsnark(t.coeff(), tmp);
        b.add_term(t.index(), tmp);
    }
}

//
// from snarklib to libsnark rank-1 constraint
//

template <typename FR,
          typename LIBSNARK_FR>
void copy_libsnark(
    const R1Constraint<FR>& a,
    libsnark::r1cs_constraint<LIBSNARK_FR>& b)
{
    libsnark::linear_combination<LIBSNARK_FR> tmpA, tmpB, tmpC;

    copy_libsnark(a.a(), tmpA);
    copy_libsnark(a.b(), tmpB);
    copy_libsnark(a.c(), tmpC);

    b = libsnark::r1cs_constraint<LIBSNARK_FR>(tmpA, tmpB, tmpC);
}

//
// from snarklib to libsnark witness vector
//

template <typename FR,
          typename LIBSNARK_FR>
void copy_libsnark(
    const R1Witness<FR>& a,
    std::vector<LIBSNARK_FR>& b,
    const std::size_t startIndex = 0,
    const std::size_t stopIndex = -1)
{
    const std::size_t L = std::min(a.size(), stopIndex);

    if (b.empty()) b.reserve(L - startIndex);

    for (std::size_t i = startIndex; i < L; ++i) {
        LIBSNARK_FR tmp;
        copy_libsnark(a[i], tmp);
        b.emplace_back(tmp);
    }
}

//
// from snarklib to libsnark rank-1 constraint system
//

template <template <typename> class SYS,
          typename FR,
          typename LIBSNARK_FR>
void copy_libsnark(
    const SYS<FR>& csA,
    const R1Witness<FR>& witnessA,
    const R1Witness<FR>& inputA,
    libsnark::r1cs_constraint_system<LIBSNARK_FR>& csB,
#ifdef USE_OLD_LIBSNARK
    libsnark::r1cs_variable_assignment<LIBSNARK_FR>& witnessB,
    libsnark::r1cs_variable_assignment<LIBSNARK_FR>& inputB)
#else
    libsnark::r1cs_auxiliary_input<LIBSNARK_FR>& witnessB,
    libsnark::r1cs_primary_input<LIBSNARK_FR>& inputB)
#endif
{
    const std::size_t numInputs = inputA.size();

#ifdef USE_OLD_LIBSNARK
    csB.num_inputs = numInputs;
    csB.num_vars = witnessA.size();
#else
    csB.primary_input_size = numInputs;
    csB.auxiliary_input_size = witnessA.size() - numInputs;
#endif

    csA.mapLambda(
        [&csB] (const R1System<FR>& a) -> bool {
            for (const auto& c : a.constraints()) {
                libsnark::r1cs_constraint<LIBSNARK_FR> tmp;
                copy_libsnark(c, tmp);
                csB.add_constraint(tmp);
            }
            return false; // do not write back to disk if HugeSystem
        });

#ifdef USE_OLD_LIBSNARK
    copy_libsnark(witnessA, witnessB);
#else
    copy_libsnark(witnessA, witnessB, numInputs);
#endif

    copy_libsnark(inputA, inputB);
}

//
// from libsnark to snarklib encoded IC (input constraint) query
//

template <typename PAIRING>
void copy_libsnark(
#ifdef USE_OLD_LIBSNARK
    const libsnark::r1cs_ppzksnark_IC_query<LIBSNARK_PPT>& a,
#else
    const libsnark::accumulation_vector<LIBSNARK_G1>& a,
#endif
    PPZK_QueryIC<PAIRING>& b)
{
    typedef typename PAIRING::G1 G1;

    G1 base;
    std::vector<G1> encoded_terms;

#ifdef USE_OLD_LIBSNARK
    copy_libsnark(a.base, base);
    copy_libsnark(a.encoded_terms, encoded_terms);
#else
    copy_libsnark(a.first, base);
    copy_libsnark(a.rest.values, encoded_terms);
#endif

    b = PPZK_QueryIC<PAIRING>(base, encoded_terms);
}

//
// from snarklib encoded IC (input constraint) query to libsnark
//

template <typename PAIRING>
void copy_libsnark(
    const PPZK_QueryIC<PAIRING>& a,
#ifdef USE_OLD_LIBSNARK
    libsnark::r1cs_ppzksnark_IC_query<LIBSNARK_PPT>& b)
#else
    libsnark::accumulation_vector<LIBSNARK_G1>& b)
#endif
{
    LIBSNARK_G1 base;
    copy_libsnark(a.base(), base);

    std::vector<LIBSNARK_G1> encoded_terms;
    copy_libsnark(a.encoded_terms(), encoded_terms);

#ifdef USE_OLD_LIBSNARK
    b = libsnark::r1cs_ppzksnark_IC_query<LIBSNARK_PPT>(base, encoded_terms);
#else
    b = libsnark::accumulation_vector<LIBSNARK_G1>(
            LIBSNARK_G1(base),                        // must be r-value
            std::vector<LIBSNARK_G1>(encoded_terms)); // must be r-value
#endif
}

//
// from libsnark to snarklib proving key
//

template <typename PAIRING>
void copy_libsnark(
    const libsnark::r1cs_ppzksnark_proving_key<LIBSNARK_PPT>& a,
    PPZK_ProvingKey<PAIRING>& b)
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    SparseVector<Pairing<G1, G1>> A_query, C_query;
    SparseVector<Pairing<G2, G1>> B_query;
    std::vector<G1> H_query, K_query;

#ifdef USE_OLD_LIBSNARK
    copy_libsnark(a.A_query, A_query);
    copy_libsnark(a.B_query, B_query);
    copy_libsnark(a.C_query, C_query);
    copy_libsnark(a.H_query, H_query);
    copy_libsnark(a.K_query, K_query);
#else
    // new and old libsnark have different query vector formats
    // new libsnark appends inhomogeneous Z values to the back
    // old libsnark prepends inhomogeneous Z values to the front
    // convert new format to old format which snarklib expects

    // A
    {
        const auto len = a.A_query.size();
        A_query.reserve(len + 2);
        G1 tmpG, tmpH;
        const auto Z = a.A_query.values.back();
        copy_libsnark(Z.g, tmpG);
        copy_libsnark(Z.h, tmpH);
        A_query.pushBack(0, Pairing<G1, G1>(tmpG, tmpH));
        A_query.pushBack(1, Pairing<G1, G1>::zero());
        A_query.pushBack(2, Pairing<G1, G1>::zero());
        copy_libsnark(a.A_query, A_query, 0, len - 1);
    }

    // B
    {
        const auto len = a.B_query.size();
        B_query.reserve(len + 2);
        G2 tmpG;
        G1 tmpH;
        const auto Z = a.B_query.values.back();
        copy_libsnark(Z.g, tmpG);
        copy_libsnark(Z.h, tmpH);
        B_query.pushBack(0, Pairing<G2, G1>::zero());
        B_query.pushBack(1, Pairing<G2, G1>(tmpG, tmpH));
        B_query.pushBack(2, Pairing<G2, G1>::zero());
        copy_libsnark(a.B_query, B_query, 0, len - 1);
    }

    // C
    {
        const auto len = a.C_query.size();
        C_query.reserve(len + 2);
        G1 tmpG, tmpH;
        const auto Z = a.C_query.values.back();
        copy_libsnark(Z.g, tmpG);
        copy_libsnark(Z.h, tmpH);
        C_query.pushBack(0, Pairing<G1, G1>::zero());
        C_query.pushBack(1, Pairing<G1, G1>::zero());
        C_query.pushBack(2, Pairing<G1, G1>(tmpG, tmpH));
        copy_libsnark(a.C_query, C_query, 0, len - 1);
    }

    // H
    copy_libsnark(a.H_query, H_query);

    // K
    {
        const auto len = a.K_query.size();
        K_query.reserve(len);
        G1 tmpG;
        copy_libsnark(a.K_query[len - 3], tmpG);
        K_query.emplace_back(tmpG);
        copy_libsnark(a.K_query[len - 2], tmpG);
        K_query.emplace_back(tmpG);
        copy_libsnark(a.K_query[len - 1], tmpG);
        K_query.emplace_back(tmpG);
        copy_libsnark(a.K_query, K_query, 0, len - 3);
    }
#endif

    b = PPZK_ProvingKey<PAIRING>(A_query,
                                 B_query,
                                 C_query,
                                 H_query,
                                 K_query);
}

//
// from snarklib to libsnark proving key
//

template <typename PAIRING>
void copy_libsnark(
    const PPZK_ProvingKey<PAIRING>& a,
    const libsnark::r1cs_constraint_system<LIBSNARK_FR>& cs,
    libsnark::r1cs_ppzksnark_proving_key<LIBSNARK_PPT>& b)
{
#ifdef USE_OLD_LIBSNARK
    copy_libsnark(a.A_query(), b.A_query);
    copy_libsnark(a.B_query(), b.B_query);
    copy_libsnark(a.C_query(), b.C_query);
    copy_libsnark(a.H_query(), b.H_query);
    copy_libsnark(a.K_query(), b.K_query);
#else
    // QAP query vectors A, B, C, K are the same size
    const std::size_t qap_num_vars = a.K_query().size() - 4;

    // A
    copy_libsnark(
        a.A_query(),
        qap_num_vars + 1,
        a.A_query().getElementForIndex(0),
        b.A_query,
        3);

    // B
    copy_libsnark(
        a.B_query(),
        qap_num_vars + 1,
        a.B_query().getElementForIndex(1),
        b.B_query,
        3);

    // C
    copy_libsnark(
        a.C_query(),
        qap_num_vars + 1,
        a.C_query().getElementForIndex(2),
        b.C_query,
        3);

    // H
    copy_libsnark(
        a.H_query(),
        b.H_query);

    // K
    copy_libsnark(a.K_query(), b.K_query, 3);
    for (std::size_t i = 0; i < 3; ++i) {
        LIBSNARK_G1 tmp;
        copy_libsnark(a.K_query()[i], tmp);
        b.K_query.emplace_back(tmp);
    }
#endif

    b.constraint_system = cs;
}

//
// from libsnark to snarklib verification key
//

template <typename PAIRING>
void copy_libsnark(
    const libsnark::r1cs_ppzksnark_verification_key<LIBSNARK_PPT>& a,
    PPZK_VerificationKey<PAIRING>& b)
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    G1 alphaB_g1, gamma_beta_g1;
    G2 alphaA_g2, alphaC_g2, gamma_g2, gamma_beta_g2, rC_Z_g2;

    copy_libsnark(a.alphaA_g2, alphaA_g2);
    copy_libsnark(a.alphaB_g1, alphaB_g1);
    copy_libsnark(a.alphaC_g2, alphaC_g2);
    copy_libsnark(a.gamma_g2, gamma_g2);
    copy_libsnark(a.gamma_beta_g1, gamma_beta_g1);
    copy_libsnark(a.gamma_beta_g2, gamma_beta_g2);
    copy_libsnark(a.rC_Z_g2, rC_Z_g2);

    PPZK_QueryIC<PAIRING> icqB;
#ifdef USE_OLD_LIBSNARK
    copy_libsnark(*a.encoded_IC_query, icqB);
#else
    copy_libsnark(a.encoded_IC_query, icqB);
#endif

    b = PPZK_VerificationKey<PAIRING>(alphaA_g2,
                                      alphaB_g1,
                                      alphaC_g2,
                                      gamma_g2,
                                      gamma_beta_g1,
                                      gamma_beta_g2,
                                      rC_Z_g2,
                                      icqB);
}

//
// from snarklib to libsnark verification key
//

template <typename PAIRING>
void copy_libsnark(
    const PPZK_VerificationKey<PAIRING>& a,
    libsnark::r1cs_ppzksnark_verification_key<LIBSNARK_PPT>& b)
{
    copy_libsnark(a.alphaA_g2(), b.alphaA_g2);
    copy_libsnark(a.alphaB_g1(), b.alphaB_g1);
    copy_libsnark(a.alphaC_g2(), b.alphaC_g2);
    copy_libsnark(a.gamma_g2(), b.gamma_g2);
    copy_libsnark(a.gamma_beta_g1(), b.gamma_beta_g1);
    copy_libsnark(a.gamma_beta_g2(), b.gamma_beta_g2);
    copy_libsnark(a.rC_Z_g2(), b.rC_Z_g2);

#ifdef USE_OLD_LIBSNARK
    b.encoded_IC_query.reset(new libsnark::r1cs_ppzksnark_IC_query<LIBSNARK_PPT>);
    copy_libsnark(a.encoded_IC_query(), *b.encoded_IC_query);
#else
    copy_libsnark(a.encoded_IC_query(), b.encoded_IC_query);
#endif
}

//
// from libsnark to snarklib proof
//

template <typename PAIRING>
void copy_libsnark(
    const libsnark::r1cs_ppzksnark_proof<LIBSNARK_PPT>& a,
    PPZK_Proof<PAIRING>& b)
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

    Pairing<G1, G1> tmpA, tmpC;
    Pairing<G2, G1> tmpB;
    G1 tmpH, tmpK;

    copy_libsnark(a.g_A, tmpA);
    copy_libsnark(a.g_B, tmpB);
    copy_libsnark(a.g_C, tmpC);
    copy_libsnark(a.g_H, tmpH);
    copy_libsnark(a.g_K, tmpK);

    b = PPZK_Proof<PAIRING>(tmpA, tmpB, tmpC, tmpH, tmpK);
}

//
// from snarklib to libsnark proof
//

template <typename PAIRING>
void copy_libsnark(
    const PPZK_Proof<PAIRING>& a,
    libsnark::r1cs_ppzksnark_proof<LIBSNARK_PPT>& b)
{
    copy_libsnark(a.A(), b.g_A);
    copy_libsnark(a.B(), b.g_B);
    copy_libsnark(a.C(), b.g_C);
    copy_libsnark(a.H(), b.g_H);
    copy_libsnark(a.K(), b.g_K);
}

} // namespace snarklib

#endif
