#ifndef _SNARKLIB_FOREIGN_LIB_HPP_
#define _SNARKLIB_FOREIGN_LIB_HPP_

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
#include /*libsnark*/ "encoding/knowledge_commitment.hpp"
#include /*libsnark*/ "r1cs/r1cs.hpp"
#else
#include /*libsnark*/ "common/default_types/ec_pp.hpp"
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

namespace snarklib {

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

template <typename UG,
          typename UH,
          typename TG,
          typename TH>
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
#ifdef USE_OLD_LIBSNARK
    const libsnark::G1_vector<libsnark::default_pp>& a,
#else
    const libsnark::G1_vector<libsnark::default_ec_pp>& a,
#endif
    std::vector<G1>& b)
{
    const std::size_t vecSize = a.size();

    if (b.empty()) b.reserve(vecSize);

    G1 tmp;

    for (std::size_t i = 0; i < vecSize; ++i) {
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
#ifdef USE_OLD_LIBSNARK
    libsnark::G1_vector<libsnark::default_pp>& b)
#else
    libsnark::G1_vector<libsnark::default_ec_pp>& b)
#endif
{
    const std::size_t vecSize = a.size();

    b.clear();
    b.reserve(vecSize);

#ifdef USE_OLD_LIBSNARK
    libsnark::G1<libsnark::default_pp> tmp;
#else
    libsnark::G1<libsnark::default_ec_pp> tmp;
#endif

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a[i], tmp);
        b.emplace_back(tmp);
    }
}

//
// from libsnark paired groups vector to snarklib sparse vector
//

template <typename G1>
void copy_libsnark(
#ifdef USE_OLD_LIBSNARK
    const libsnark::G1G1_knowledge_commitment_vector<libsnark::default_pp>& a,
#else
    const libsnark::knowledge_commitment_vector<libsnark::G1<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>& a,
#endif
    SparseVector<Pairing<G1, G1>>& b)
{
    const std::size_t
        vecSize = a.values.size(),
        offset = b.size();

    if (b.empty()) b.reserve(vecSize);

    G1 tmpG, tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].g, tmpG);
        copy_libsnark(a.values[i].h, tmpH);

        b.pushBack(
            a.indices[i] + offset,
            Pairing<G1, G1>(tmpG, tmpH));
    }
}

template <typename G2,
          typename G1>
void copy_libsnark(
#ifdef USE_OLD_LIBSNARK
    const libsnark::G2G1_knowledge_commitment_vector<libsnark::default_pp>& a,
#else
    const libsnark::knowledge_commitment_vector<libsnark::G2<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>& a,
#endif
    SparseVector<Pairing<G2, G1>>& b)
{
    const std::size_t
        vecSize = a.values.size(),
        offset = b.size();

    if (b.empty()) b.reserve(vecSize);

    G2 tmpG;
    G1 tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].g, tmpG);
        copy_libsnark(a.values[i].h, tmpH);

        b.pushBack(
            a.indices[i] + offset,
            Pairing<G2, G1>(tmpG, tmpH));
    }
}

//
// from snarklib sparse vector to libsnark paired groups vector
//

#ifdef USE_OLD_LIBSNARK
template <typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G1, G1>>& a,
    libsnark::G1G1_knowledge_commitment_vector<libsnark::default_pp>& b)
{
    const std::size_t vecSize = a.values.size();

    b.values.clear();
    b.values.reserve(vecSize);
    b.indices.clear();
    b.indices.reserve(vecSize);
    b.is_sparse = true;
    b.original_size = vecSize;

    libsnark::G1<libsnark::default_pp> tmpG, tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].G(), tmpG);
        copy_libsnark(a.values[i].H(), tmpH);

        b.values.emplace_back(
            libsnark::knowledge_commitment<libsnark::G1<libsnark::default_pp>, libsnark::G1<libsnark::default_pp>>(tmpG, tmpH));

        b.indices.push_back(a.indices[i]);
    }
}

template <typename G2,
          typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G2, G1>>& a,
    libsnark::G2G1_knowledge_commitment_vector<libsnark::default_pp>& b)
{
    const std::size_t vecSize = a.values.size();

    b.values.clear();
    b.values.reserve(vecSize);
    b.indices.clear();
    b.indices.reserve(vecSize);
    b.is_sparse = true;
    b.original_size = vecSize;

    libsnark::G2<libsnark::default_pp> tmpG;
    libsnark::G1<libsnark::default_pp> tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].G(), tmpG);
        copy_libsnark(a.values[i].H(), tmpH);

        b.values.emplace_back(
            libsnark::knowledge_commitment<libsnark::G2<libsnark::default_pp>, libsnark::G1<libsnark::default_pp>>(tmpG, tmpH));

        b.indices.push_back(a.indices[i]);
    }
}
#else
template <typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G1, G1>>& a,
    libsnark::knowledge_commitment_vector<libsnark::G1<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>& b)
{
    const std::size_t vecSize = a.values.size();

    b.indices.clear();
    b.indices.reserve(vecSize);
    b.values.clear();
    b.values.reserve(vecSize);
    b.domain_size_ = vecSize;

    libsnark::G1<libsnark::default_ec_pp> tmpG, tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].G(), tmpG);
        copy_libsnark(a.values[i].H(), tmpH);

        b.indices.push_back(a.indices[i]);

        b.values.emplace_back(
            libsnark::knowledge_commitment<libsnark::G1<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>(tmpG, tmpH));
    }
}

template <typename G2,
          typename G1>
void copy_libsnark(
    const SparseVector<Pairing<G2, G1>>& a,
    libsnark::knowledge_commitment_vector<libsnark::G2<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>& b)
{
    const std::size_t vecSize = a.values.size();

    b.indices.clear();
    b.indices.reserve(vecSize);
    b.values.clear();
    b.values.reserve(vecSize);
    b.domain_size_ = vecSize;

    libsnark::G2<libsnark::default_ec_pp> tmpG;
    libsnark::G1<libsnark::default_ec_pp> tmpH;

    for (std::size_t i = 0; i < vecSize; ++i) {
        copy_libsnark(a.values[i].G(), tmpG);
        copy_libsnark(a.values[i].H(), tmpH);

        b.indices.push_back(a.indices[i]);

        b.values.emplace_back(
            libsnark::knowledge_commitment<libsnark::G2<libsnark::default_ec_pp>, libsnark::G1<libsnark::default_ec_pp>>(tmpG, tmpH));
    }
}
#endif

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
        });

#ifdef USE_OLD_LIBSNARK
    for (const auto& w : *witnessA) {
        LIBSNARK_FR tmp;
        copy_libsnark(w, tmp);
        witnessB.push_back(tmp);
    }
#else
    for (std::size_t i = numInputs; i < witnessA.size(); ++i) {
        LIBSNARK_FR tmp;
        copy_libsnark(witnessA[i], tmp);
        witnessB.push_back(tmp);
    }
#endif

    for (const auto& w : *inputA) {
        LIBSNARK_FR tmp;
        copy_libsnark(w, tmp);
        inputB.push_back(tmp);
    }
}

} // namespace snarklib

#endif
