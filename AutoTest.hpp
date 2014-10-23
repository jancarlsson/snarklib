#ifndef _SNARKLIB_AUTO_TEST_HPP_
#define _SNARKLIB_AUTO_TEST_HPP_

#include <cassert>
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <memory>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include "algebra/curves/alt_bn128/alt_bn128_g1.hpp"
#include "algebra/curves/alt_bn128/alt_bn128_g2.hpp"
#include "algebra/curves/edwards/edwards_g1.hpp"
#include "algebra/curves/edwards/edwards_g2.hpp"
#include "algebra/fields/bigint.hpp"
#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp2.hpp"
#include "algebra/fields/fp3.hpp"
#include "algebra/fields/fp6_2over3.hpp"
#include "algebra/fields/fp6_3over2.hpp"
#include "algebra/fields/fp12_2over3over2.hpp"
#include "BigInt.hpp"
#include "common/types.hpp"
#include "EC.hpp"
#include "encoding/knowledge_commitment.hpp"
#include "Field.hpp"
#include "FpModel.hpp"
#include "Group.hpp"
#include "Pairing.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// autotest base class
//

class AutoTest
{
public:
    ~AutoTest() = default;

    std::string testName() const {
        std::stringstream ss;
        ss << typeid(*this).name() << m_testNameSuffix;
        return ss.str();
    }

    std::size_t testNumber() const {
        return m_testCounter;
    }

    virtual void runTest() = 0;

    bool testPass() const {
        return m_testPass;
    }

    void testLog(std::ostream& out) const {
        out << testNumber() << "\t"
            << (testPass() ? "PASS" : "FAIL") << "\t"
            << testName() << std::endl;
    }

protected:
    AutoTest()
        : m_testPass(true)
    {
        static std::size_t testCounter = 0;
        m_testCounter = testCounter++;
    }

    template <typename A>
    AutoTest(const A& a)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B>
    AutoTest(const A& a, const B& b)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B, typename C>
    AutoTest(const A& a, const B& b, const C& c)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b << " " << c;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B, typename C, typename D>
    AutoTest(const A& a, const B& b, const C& c, const D& d)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b << " " << c << " " << d;
        m_testNameSuffix = ss.str();
    }

    bool checkPass(const bool v) {
        m_testPass &= v;

        return v;
    }

    // bigint == BigInt
    template <mp_size_t N>
    bool sameData(const libsnark::bigint<N>& a,
                  const BigInt<N>& b) const
    {
        for (std::size_t i = 0; i < N; ++i) {
            if (a.data[i] != b.data()[i])
                return false;
        }

        return true;
    }

    // Fp_model == FpModel
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp_model<N, MODULUS_A>& a,
                  const FpModel<N, MODULUS_B>& b) const
    {
        return sameData(a.as_bigint(), b.asBigInt());
    }

    // Fp_model == Field<FpModel>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp_model<N, MODULUS_A>& a,
                  const Field<FpModel<N, MODULUS_B>>& b) const
    {
        return sameData(a, b[0]);
    }

    // Fp2_model == Field<FpModel, 2>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp2_model<N, MODULUS_A>& a,
                  const Field<FpModel<N, MODULUS_B>, 2>& b) const
    {
        return
            sameData(a.c0, b[0]) &&
            sameData(a.c1, b[1]);
    }

    // Fp3_model == Field<FpModel, 3>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp3_model<N, MODULUS_A>& a,
                  const Field<FpModel<N, MODULUS_B>, 3>& b) const
    {
        return
            sameData(a.c0, b[0]) &&
            sameData(a.c1, b[1]) &&
            sameData(a.c2, b[2]);
    }

    // Fp6_2over3_model == Field<Field<FpModel, 3>, 2>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp6_2over3_model<N, MODULUS_A>& a,
                  const Field<Field<FpModel<N, MODULUS_B>, 3>, 2>& b) const
    {
        return
            sameData(a.c0, b[0]) &&
            sameData(a.c1, b[1]);
    }

    // Fp6_3over2_model == Field<Field<FpModel, 2>, 3>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp6_3over2_model<N, MODULUS_A>& a,
                  const Field<Field<FpModel<N, MODULUS_B>, 2>, 3>& b) const
    {
        return
            sameData(a.c0, b[0]) &&
            sameData(a.c1, b[1]) &&
            sameData(a.c2, b[2]);
    }

    // Fp12_2over3over2_model == Field<Field<Field<FpModel, 2>, 3>, 2>
    template <mp_size_t N, const libsnark::bigint<N>& MODULUS_A, const BigInt<N>& MODULUS_B>
    bool sameData(const libsnark::Fp12_2over3over2_model<N, MODULUS_A>& a,
                  const Field<Field<Field<FpModel<N, MODULUS_B>, 2>, 3>, 2>& b) const
    {
        return
            sameData(a.c0, b[0]) &&
            sameData(a.c1, b[1]);
    }

    // alt_bn128_G1 == BN128::Groups<>::G1
    template <typename GROUP>
    bool sameData(const libsnark::alt_bn128_G1& a, const GROUP& b) const {
        return
            sameData(a.X, b.x()) &&
            sameData(a.Y, b.y()) &&
            sameData(a.Z, b.z());
    }

    // alt_bn128_G2 == BN128::Groups<>::G2
    template <typename GROUP>
    bool sameData(const libsnark::alt_bn128_G2& a, const GROUP& b) const {
        return
            sameData(a.X, b.x()) &&
            sameData(a.Y, b.y()) &&
            sameData(a.Z, b.z());
    }

    // edwards_G1 == Edwards::Groups<>::G1
    template <typename GROUP>
    bool sameData(const libsnark::edwards_G1& a, const GROUP& b) const {
        return
            sameData(a.X, b.x()) &&
            sameData(a.Y, b.y()) &&
            sameData(a.Z, b.z());
    }

    // edwards_G2 == Edwards::Groups<>::G2
    template <typename GROUP>
    bool sameData(const libsnark::edwards_G2& a, const GROUP& b) const {
        return
            sameData(a.X, b.x()) &&
            sameData(a.Y, b.y()) &&
            sameData(a.Z, b.z());
    }

    // knowledge_commitment<> == Pairing<>
    template <typename UG, typename UH, typename TG, typename TH>
    bool sameData(const libsnark::knowledge_commitment<UG, UH>& a,
                  const Pairing<TG, TH>& b) const {
        return
            sameData(a.g, b.G()) &&
            sameData(a.h, b.H());
    }

    // vector<> == vector<>
    template <typename T, typename U>
    bool sameData(const std::vector<T>& a,
                  const std::vector<U>& b) const {
        if (a.size() != b.size()) {
            return false;
        }

        for (std::size_t i = 0; i < a.size(); ++i) {
            if (! sameData(a[i], b[i]))
                return false;
        }

        return true;
    }

    // BigInt -> bigint
    template <mp_size_t N>
    void copyData(const BigInt<N>& a, libsnark::bigint<N>& b) {
        std::stringstream ssA;
        ssA << a;
        std::stringstream ssB(ssA.str());
        ssB >> b;
    }

    // bigint -> BigInt
    template <mp_size_t N>
    void copyData(const libsnark::bigint<N>& a, BigInt<N>& b) {
        std::stringstream ssA;
        ssA << a;
        std::stringstream ssB(ssA.str());
        ssB >> b;
    }

#define DEFN_COPY_DATA_F2M(A, B)                        \
    template <mp_size_t N,                              \
              const BigInt<N>& MA,                      \
              const libsnark::bigint<N>& MB>            \
    void copyData(const A & a, libsnark:: B & b) {      \
        std::stringstream ssA;                          \
        ssA << a;                                       \
        std::stringstream ssB(ssA.str());               \
        ssB >> b;                                       \
    }

#define DEFN_COPY_DATA_M2F(A, B)                        \
    template <mp_size_t N,                              \
              const libsnark::bigint<N>& MA,            \
              const BigInt<N>& MB>                      \
    void copyData(const libsnark:: A & a, B & b) {      \
        std::stringstream ssA;                          \
        ssA << a;                                       \
        std::stringstream ssB(ssA.str());               \
        ssB >> b;                                       \
    }

#define DEFN_COPY_DATA_G2M(B)                           \
    template <typename GROUP>                           \
    void copyData(const GROUP& a, libsnark:: B & b) {   \
        copyData(a.x(), b.X);                           \
        copyData(a.y(), b.Y);                           \
        copyData(a.z(), b.Z);                           \
    }

#define DEFN_COPY_DATA_M2G(A)                           \
    template <typename GROUP>                           \
    void copyData(const libsnark:: A & a, GROUP& b) {   \
    typename GROUP::BaseField tmp;                      \
        copyData(a.X, tmp);                             \
        b.x(tmp);                                       \
        copyData(a.Y, tmp);                             \
        b.y(tmp);                                       \
        copyData(a.Z, tmp);                             \
        b.z(tmp);                                       \
    }

#define COMMA ,

    // from snarklib Field<> to libsnark FpX_model
    DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA>>, Fp_model<N COMMA MB>)
    DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA> COMMA 2>, Fp2_model<N COMMA MB>)
    DEFN_COPY_DATA_F2M(Field<FpModel<N COMMA MA> COMMA 3>, Fp3_model<N COMMA MB>)
    //DEFN_COPY_DATA_F2M(Field<Field<FpModel<N COMMA MA> COMMA 3> COMMA 2>, Fp6_2over3_model<N COMMA MB>)
    DEFN_COPY_DATA_F2M(Field<Field<FpModel<N COMMA MA> COMMA 2> COMMA 3>, Fp6_3over2_model<N COMMA MB>)
    DEFN_COPY_DATA_F2M(Field<Field<Field<FpModel<N COMMA MA> COMMA 2> COMMA 3> COMMA 2>, Fp12_2over3over2_model<N COMMA MB>)

    // from libsnark FpX_model to snarklib Field<>
    DEFN_COPY_DATA_M2F(Fp_model<N COMMA MA>, Field<FpModel<N COMMA MB>>)
    DEFN_COPY_DATA_M2F(Fp2_model<N COMMA MA>, Field<FpModel<N COMMA MB> COMMA 2>)
    DEFN_COPY_DATA_M2F(Fp3_model<N COMMA MA>, Field<FpModel<N COMMA MB> COMMA 3>)
    //DEFN_COPY_DATA_M2F(Fp6_2over3_model<N COMMA MA>, Field<Field<FpModel<N COMMA MB> COMMA 3> COMMA 2>)
    DEFN_COPY_DATA_M2F(Fp6_3over2_model<N COMMA MA>, Field<Field<FpModel<N COMMA MB> COMMA 2> COMMA 3>)
    DEFN_COPY_DATA_M2F(Fp12_2over3over2_model<N COMMA MA>, Field<Field<Field<FpModel<N COMMA MB> COMMA 2> COMMA 3> COMMA 2>)

    // special case for snarklib Field<> to libsnark Fp32 model
    template <mp_size_t N, const BigInt<N>& MA, const libsnark::bigint<N>& MB>
    void copyData(const Field<Field<FpModel<N, MA>, 3>, 2>& a,
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

    // special case for libsnark Fp32 model to snarklib Field<>
    template <mp_size_t N, const libsnark::bigint<N>& MA, const BigInt<N>& MB>
    void copyData(const libsnark::Fp6_2over3_model<N, MA>& a,
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

    // from snarklib Group<> to libsnark elliptic curve group
    DEFN_COPY_DATA_G2M(alt_bn128_G1)
    DEFN_COPY_DATA_G2M(alt_bn128_G2)
    DEFN_COPY_DATA_G2M(edwards_G1)
    DEFN_COPY_DATA_G2M(edwards_G2)

    // from libsnark elliptic curve group to snarklib Group<>
    DEFN_COPY_DATA_M2G(alt_bn128_G1)
    DEFN_COPY_DATA_M2G(alt_bn128_G2)
    DEFN_COPY_DATA_M2G(edwards_G1)
    DEFN_COPY_DATA_M2G(edwards_G2)

#undef DEFN_COPY_DATA_F2M
#undef DEFN_COPY_DATA_M2F
#undef DEFN_COPY_DATA_G2M
#undef DEFN_COPY_DATA_M2G
#undef COMMA

    // from libsnark group vector to snarklib vector
    template <typename G1>
    void copyData(
        const libsnark::G1_vector<libsnark::default_pp>& a,
        std::vector<G1>& b)
    {
        const std::size_t vecSize = a.size();

        b.clear();
        b.reserve(vecSize);

        G1 tmp;

        for (std::size_t i = 0; i < vecSize; ++i) {
            copyData(a[i], tmp);
            b.emplace_back(tmp);
        }
    }

    // from libsnark paired groups vector to snarklib sparse vector
    template <typename G1>
    void copyData(
        const libsnark::G1G1_knowledge_commitment_vector<libsnark::default_pp>& a,
        SparseVector<Pairing<G1, G1>>& b)
    {
        const std::size_t vecSize = a.values.size();

        b.clear();
        b.reserve(vecSize);

        G1 tmpG, tmpH;

        for (std::size_t i = 0; i < vecSize; ++i) {
            copyData(a.values[i].g, tmpG);
            copyData(a.values[i].h, tmpH);

            b.pushBack(
                a.indices[i],
                Pairing<G1, G1>(tmpG, tmpH));
        }
    }

    // from libsnark paired groups vector to snarklib sparse vector
    template <typename G2, typename G1>
    void copyData(
        const libsnark::G2G1_knowledge_commitment_vector<libsnark::default_pp>& a,
        SparseVector<Pairing<G2, G1>>& b)
    {
        const std::size_t vecSize = a.values.size();

        b.clear();
        b.reserve(vecSize);

        G2 tmpG;
        G1 tmpH;

        for (std::size_t i = 0; i < vecSize; ++i) {
            copyData(a.values[i].g, tmpG);
            copyData(a.values[i].h, tmpH);

            b.pushBack(
                a.indices[i],
                Pairing<G2, G1>(tmpG, tmpH));
        }
    }

    // vector<>
    template <typename T>
    void printData(std::ostream& out,
                   const std::string& prefix,
                   const std::vector<T>& v) const
    {
        for (std::size_t i = 0; i < v.size(); ++i) {
            out << prefix << "[" << i << "] = " << v[i] << std::endl;
        }
    }

private:
    std::string m_testNameSuffix;
    std::size_t m_testCounter;
    bool m_testPass;
};

////////////////////////////////////////////////////////////////////////////////
// battery of autotests
//

class AutoTestBattery
{
public:
    AutoTestBattery() = default;

    std::size_t testCount() const {
        return m_testVector.size();
    }

    void addTest(AutoTest* ptr) {
        m_testVector.push_back(std::unique_ptr<AutoTest>(ptr));
    }

    bool runTest() {
        bool allPass = true;

        for (const auto& a : m_testVector) {
            a->runTest();
            allPass &= a->testPass();
        }

        return allPass;
    }

    void testLog(std::ostream& out) const {
        for (const auto& a : m_testVector) {
            a->testLog(out);
        }
    }

private:
    std::vector<std::unique_ptr<AutoTest>> m_testVector;
};

////////////////////////////////////////////////////////////////////////////////
// convenient functions
//

std::string randomBase10(std::random_device& rd, const mp_size_t N) {
    std::stringstream ss;
    for (size_t i = 0; i < N; ++i) {
        ss << rd();
    }

    return ss.str();
}

std::string uniformBase10(const unsigned long low, const unsigned long high) {
    assert(low <= high);

    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned long> distribution(low, high);

    std::stringstream ss;
    ss << distribution(generator);

    return ss.str();
}

std::string sparseUniformBase10(const unsigned long low, const unsigned long high) {
    assert(low <= high);

    std::default_random_engine generator;

    std::uniform_int_distribution<unsigned long>
        selector(0, 2),
        distribution(low, high);

    std::stringstream ss;

    switch (selector(generator)) {
    case (0) :
        ss << 0;
        break;

    case (1) :
        ss << 1;
        break;

    default:
        ss << distribution(generator);
        break;
    }

    return ss.str();
}

template <mp_size_t N>
libsnark::bigint<N> to_bigint(const std::string& base10) {
    return libsnark::bigint<N>(base10.c_str());
}

} // namespace snarklib

#endif
