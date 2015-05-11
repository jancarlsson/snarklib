#ifndef _SNARKLIB_AUTOTEST_MARSHALLING_HPP_
#define _SNARKLIB_AUTOTEST_MARSHALLING_HPP_

#include <cstdint>
#include <gmp.h>
#include <sstream>

#include "snarklib/AutoTest.hpp"
#include "snarklib/AuxSTL.hpp"
#include "snarklib/BigInt.hpp"
#include "snarklib/Field.hpp"
#include "snarklib/Group.hpp"
#include "snarklib/Pairing.hpp"
#include "snarklib/PPZK_keypair.hpp"
#include "snarklib/PPZK_keystruct.hpp"
#include "snarklib/PPZK_proof.hpp"
#include "snarklib/PPZK_query.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// roundtrip for SparseVector<Pairing<GA, GB>>
//

template <typename GA, typename GB>
class AutoTest_Marshal_SparseVectorPairing : public AutoTest
{
public:
    AutoTest_Marshal_SparseVectorPairing(const std::size_t numberElems,
                                         const std::size_t startIndex)
        : AutoTest(numberElems, startIndex)
    {
        randomSparseVector(m_A, numberElems, startIndex);
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(oss);

        SparseVector<Pairing<GA, GB>> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(m_A == B);
    }

private:
    SparseVector<Pairing<GA, GB>> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for SparseVector<Pairing<GA, GB>> using raw data
//

template <typename GA, typename GB>
class AutoTest_Marshal_SparseVectorPairing_raw : public AutoTest
{
public:
    AutoTest_Marshal_SparseVectorPairing_raw(const std::size_t numberElems,
                                             const std::size_t startIndex)
        : AutoTest(numberElems, startIndex)
    {
        randomSparseVector(m_A, numberElems, startIndex);
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(
            oss,
            [] (std::ostream& o, const Pairing<GA, GB>& a) {
                a.marshal_out_raw(o);
            });

        SparseVector<Pairing<GA, GB>> B;

        std::stringstream iss(oss.str());
        checkPass(
            B.marshal_in(
                iss,
                [] (std::istream& i, Pairing<GA, GB>& a) {
                    return a.marshal_in_raw(i);
                }));

        checkPass(m_A == B);
    }

private:
    SparseVector<Pairing<GA, GB>> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for SparseVector<Pairing<GA, GB>> using affine representation
//

template <typename GA, typename GB>
class AutoTest_Marshal_SparseVectorPairing_special : public AutoTest
{
public:
    AutoTest_Marshal_SparseVectorPairing_special(const std::size_t numberElems,
                                                 const std::size_t startIndex)
        : AutoTest(numberElems, startIndex)
    {
        randomSparseVector(m_A, numberElems, startIndex);
        batchSpecial(m_A);
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(
            oss,
            [] (std::ostream& o, const Pairing<GA, GB>& a) {
                a.marshal_out_special(o);
            });

        SparseVector<Pairing<GA, GB>> B;

        std::stringstream iss(oss.str());
        checkPass(
            B.marshal_in(
                iss,
                [] (std::istream& i, Pairing<GA, GB>& a) {
                    return a.marshal_in_special(i);
                }));

        checkPass(m_A == B);
    }

private:
    SparseVector<Pairing<GA, GB>> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for SparseVector<Pairing<GA, GB>> using raw affine representation
//

template <typename GA, typename GB>
class AutoTest_Marshal_SparseVectorPairing_rawspecial : public AutoTest
{
public:
    AutoTest_Marshal_SparseVectorPairing_rawspecial(const std::size_t numberElems,
                                                    const std::size_t startIndex)
        : AutoTest(numberElems, startIndex)
    {
        randomSparseVector(m_A, numberElems, startIndex);
        batchSpecial(m_A);
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(
            oss,
            [] (std::ostream& o, const Pairing<GA, GB>& a) {
                a.marshal_out_rawspecial(o);
            });

        SparseVector<Pairing<GA, GB>> B;

        std::stringstream iss(oss.str());
        checkPass(
            B.marshal_in(
                iss,
                [] (std::istream& i, Pairing<GA, GB>& a) {
                    return a.marshal_in_rawspecial(i);
                }));

        checkPass(m_A == B);
    }

private:
    SparseVector<Pairing<GA, GB>> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for BigInt<>, Field<>, Group<>
//

template <typename T>
class AutoTest_Marshal_BFG : public AutoTest
{
public:
    AutoTest_Marshal_BFG()
        : AutoTest(),
          m_A(T::random())
    {}

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(oss);

        T B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(m_A == B);
    }

private:
    const T m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for BigInt<>, Field<>, Group<> using raw data
//

template <typename T>
class AutoTest_Marshal_BFG_raw : public AutoTest
{
public:
    AutoTest_Marshal_BFG_raw()
        : AutoTest(),
          m_A(T::random())
    {}

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out_raw(oss);

        T B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_raw(iss));

        checkPass(m_A == B);
    }

private:
    const T m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for vectors of Group<>
//

template <typename T>
class AutoTest_Marshal_vectorG : public AutoTest
{
public:
    AutoTest_Marshal_vectorG(const std::size_t N)
        : AutoTest(),
          m_A(N)
    {
        for (auto& a : m_A)
            a = T::random();
    }

    void runTest() {
        std::stringstream oss;
        marshal_out(oss, m_A);

        std::vector<T> B;

        std::stringstream iss(oss.str());
        checkPass(marshal_in(iss, B));

        checkPass(m_A == B);
    }

private:
    std::vector<T> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for vectors of Group<> using raw data
//

template <typename T>
class AutoTest_Marshal_vectorG_raw : public AutoTest
{
public:
    AutoTest_Marshal_vectorG_raw(const std::size_t N)
        : AutoTest(),
          m_A(N)
    {
        for (auto& a : m_A)
            a = T::random();
    }

    void runTest() {
        std::stringstream oss;
        marshal_out_raw(oss, m_A);

        std::vector<T> B;

        std::stringstream iss(oss.str());
        checkPass(marshal_in_raw(iss, B));

        checkPass(m_A == B);
    }

private:
    std::vector<T> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for Pairing<GA, GB>
//

template <typename GA, typename GB>
class AutoTest_Marshal_Pairing : public AutoTest
{
public:
    AutoTest_Marshal_Pairing()
        : AutoTest(),
          m_A(GA::random(), GB::random())
    {}

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(oss);

        Pairing<GA, GB> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(m_A == B);
    }

private:
    const Pairing<GA, GB> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for Pairing<GA, GB> using raw data
//

template <typename GA, typename GB>
class AutoTest_Marshal_Pairing_raw : public AutoTest
{
public:
    AutoTest_Marshal_Pairing_raw()
        : AutoTest(),
          m_A(GA::random(), GB::random())
    {}

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out_raw(oss);

        Pairing<GA, GB> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_raw(iss));

        checkPass(m_A == B);
    }

private:
    const Pairing<GA, GB> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for Pairing<GA, GB> using affine representation
//

template <typename GA, typename GB>
class AutoTest_Marshal_Pairing_special : public AutoTest
{
public:
    AutoTest_Marshal_Pairing_special()
        : AutoTest(),
          m_A(GA::random(), GB::random())
    {
        m_A.toSpecial();
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out_special(oss);

        Pairing<GA, GB> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_special(iss));

        checkPass(m_A == B);
    }

private:
    Pairing<GA, GB> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for Pairing<GA, GB> using raw affine representation
//

template <typename GA, typename GB>
class AutoTest_Marshal_Pairing_rawspecial : public AutoTest
{
public:
    AutoTest_Marshal_Pairing_rawspecial()
        : AutoTest(),
          m_A(GA::random(), GB::random())
    {
        m_A.toSpecial();
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out_rawspecial(oss);

        Pairing<GA, GB> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_rawspecial(iss));

        checkPass(m_A == B);
    }

private:
    Pairing<GA, GB> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_ProvingKey<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_ProvingKey : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_ProvingKey(const std::size_t numberElems,
                                const std::size_t startIndex)
        : AutoTest(numberElems, startIndex),
          m_numberElems(numberElems),
          m_startIndex(startIndex)
    {}

    void runTest() {
        SparseVector<Pairing<G1, G1>> A_query;
        SparseVector<Pairing<G2, G1>> B_query;
        SparseVector<Pairing<G1, G1>> C_query;
        std::vector<G1> H_query;
        std::vector<G1> K_query;

        randomSparseVector(A_query, m_numberElems, m_startIndex);
        randomSparseVector(B_query, m_numberElems, m_startIndex);
        randomSparseVector(C_query, m_numberElems, m_startIndex);
        randomVector(H_query, m_numberElems);
        randomVector(K_query, m_numberElems);

        const PPZK_ProvingKey<PAIRING> A(A_query,
                                         B_query,
                                         C_query,
                                         H_query,
                                         K_query);

        std::stringstream oss;
        A.marshal_out(oss);

        PPZK_ProvingKey<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems, m_startIndex;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_ProvingKey<PAIRING> using raw data
//

template <typename PAIRING>
class AutoTest_Marshal_ProvingKey_raw : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_ProvingKey_raw(const std::size_t numberElems,
                                    const std::size_t startIndex)
        : AutoTest(numberElems, startIndex),
          m_numberElems(numberElems),
          m_startIndex(startIndex)
    {}

    void runTest() {
        SparseVector<Pairing<G1, G1>> A_query;
        SparseVector<Pairing<G2, G1>> B_query;
        SparseVector<Pairing<G1, G1>> C_query;
        std::vector<G1> H_query;
        std::vector<G1> K_query;

        randomSparseVector(A_query, m_numberElems, m_startIndex);
        randomSparseVector(B_query, m_numberElems, m_startIndex);
        randomSparseVector(C_query, m_numberElems, m_startIndex);
        randomVector(H_query, m_numberElems);
        randomVector(K_query, m_numberElems);

        const PPZK_ProvingKey<PAIRING> A(A_query,
                                         B_query,
                                         C_query,
                                         H_query,
                                         K_query);

        std::stringstream oss;
        A.marshal_out_raw(oss);

        PPZK_ProvingKey<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_raw(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems, m_startIndex;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_QueryIC<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_QueryIC : public AutoTest
{
    typedef typename PAIRING::G1 G1;

public:
    AutoTest_Marshal_QueryIC(const std::size_t numberElems)
        : AutoTest(numberElems),
          m_numberElems(numberElems)
    {}

    void runTest() {
        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_QueryIC<PAIRING> A(G1::random(),
                                      encoded_terms);

        std::stringstream oss;
        A.marshal_out(oss);

        PPZK_QueryIC<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_QueryIC<PAIRING> using raw data
//

template <typename PAIRING>
class AutoTest_Marshal_QueryIC_raw : public AutoTest
{
    typedef typename PAIRING::G1 G1;

public:
    AutoTest_Marshal_QueryIC_raw(const std::size_t numberElems)
        : AutoTest(numberElems),
          m_numberElems(numberElems)
    {}

    void runTest() {
        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_QueryIC<PAIRING> A(G1::random(),
                                      encoded_terms);

        std::stringstream oss;
        A.marshal_out_raw(oss);

        PPZK_QueryIC<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_raw(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_VerificationKey<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_VerificationKey : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_VerificationKey(const std::size_t numberElems)
        : AutoTest(numberElems),
          m_numberElems(numberElems)
    {}

    void runTest() {
        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_VerificationKey<PAIRING>
            A(G2::random(),
              G1::random(),
              G2::random(),
              G2::random(),
              G1::random(),
              G2::random(),
              G2::random(),
              PPZK_QueryIC<PAIRING>(G1::random(), encoded_terms));

        std::stringstream oss;
        A.marshal_out(oss);

        PPZK_VerificationKey<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_VerificationKey<PAIRING> using raw data
//

template <typename PAIRING>
class AutoTest_Marshal_VerificationKey_raw : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_VerificationKey_raw(const std::size_t numberElems)
        : AutoTest(numberElems),
          m_numberElems(numberElems)
    {}

    void runTest() {
        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_VerificationKey<PAIRING>
            A(G2::random(),
              G1::random(),
              G2::random(),
              G2::random(),
              G1::random(),
              G2::random(),
              G2::random(),
              PPZK_QueryIC<PAIRING>(G1::random(), encoded_terms));

        std::stringstream oss;
        A.marshal_out_raw(oss);

        PPZK_VerificationKey<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in_raw(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_Keypair<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_Keypair : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_Keypair(const std::size_t numberElems,
                             const std::size_t startIndex)
        : AutoTest(numberElems, startIndex),
          m_numberElems(numberElems),
          m_startIndex(startIndex)
    {}

    void runTest() {
        SparseVector<Pairing<G1, G1>> A_query;
        SparseVector<Pairing<G2, G1>> B_query;
        SparseVector<Pairing<G1, G1>> C_query;
        std::vector<G1> H_query;
        std::vector<G1> K_query;

        randomSparseVector(A_query, m_numberElems, m_startIndex);
        randomSparseVector(B_query, m_numberElems, m_startIndex);
        randomSparseVector(C_query, m_numberElems, m_startIndex);
        randomVector(H_query, m_numberElems);
        randomVector(K_query, m_numberElems);

        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_Keypair<PAIRING>
            A(PPZK_ProvingKey<PAIRING>(
                  A_query,
                  B_query,
                  C_query,
                  H_query,
                  K_query),
              PPZK_VerificationKey<PAIRING>(
                  G2::random(),
                  G1::random(),
                  G2::random(),
                  G2::random(),
                  G1::random(),
                  G2::random(),
                  G2::random(),
                  PPZK_QueryIC<PAIRING>(G1::random(), encoded_terms)));

        std::stringstream oss;
        A.marshal_out(oss);

        PPZK_Keypair<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(A == B);
    }

private:
    const std::size_t m_numberElems, m_startIndex;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for PPZK_Proof<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_Proof : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_Marshal_Proof()
        : AutoTest(),
          m_A(Pairing<G1, G1>(G1::random(), G1::random()),
              Pairing<G2, G1>(G2::random(), G1::random()),
              Pairing<G1, G1>(G1::random(), G1::random()),
              G1::random(),
              G1::random())
    {}

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(oss);

        PPZK_Proof<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(m_A == B);
    }

private:
    const PPZK_Proof<PAIRING> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// roundtrip for R1Witness<FIELD>
//

template <typename FIELD>
class AutoTest_Marshal_R1Witness : public AutoTest
{
public:
    AutoTest_Marshal_R1Witness(const std::size_t numberElems)
        : AutoTest(numberElems)
    {
        for (std::size_t i = 0; i < numberElems; ++i) {
            m_A.assignVar(R1Variable<FIELD>(i + 1), FIELD::random());
        }
    }

    void runTest() {
        std::stringstream oss;
        m_A.marshal_out(oss);

        R1Witness<FIELD> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));

        checkPass(m_A == B);
    }

private:
    R1Witness<FIELD> m_A;
};

} // namespace snarklib

#endif
