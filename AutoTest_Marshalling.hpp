#ifndef _SNARKLIB_AUTOTEST_MARSHALLING_HPP_
#define _SNARKLIB_AUTOTEST_MARSHALLING_HPP_

#include <cstdint>
#include <gmp.h>
#include <sstream>
#include "AutoTest.hpp"
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "Field.hpp"
#include "Group.hpp"
#include "Pairing.hpp"
#include "PPZK.hpp"

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
// roundtrip for PPZK_IC_Query<PAIRING>
//

template <typename PAIRING>
class AutoTest_Marshal_IC_Query : public AutoTest
{
    typedef typename PAIRING::G1 G1;

public:
    AutoTest_Marshal_IC_Query(const std::size_t numberElems)
        : AutoTest(numberElems),
          m_numberElems(numberElems)
    {}

    void runTest() {
        std::vector<G1> encoded_terms;
        randomVector(encoded_terms, m_numberElems);

        const PPZK_IC_Query<PAIRING> A(G1::random(),
                                       encoded_terms);

        std::stringstream oss;
        A.marshal_out(oss);

        PPZK_IC_Query<PAIRING> B;

        std::stringstream iss(oss.str());
        checkPass(B.marshal_in(iss));
        
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
              PPZK_IC_Query<PAIRING>(G1::random(), encoded_terms));

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
                  PPZK_IC_Query<PAIRING>(G1::random(), encoded_terms)));

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
