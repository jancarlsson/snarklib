#ifndef _SNARKLIB_PPZK_KEYPAIR_HPP_
#define _SNARKLIB_PPZK_KEYPAIR_HPP_

#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include "AuxSTL.hpp"
#include "Group.hpp"
#include "Pairing.hpp"
#include "PPZK_keystruct.hpp"
#include "PPZK_query.hpp"
#include "PPZK_randomness.hpp"
#include "ProgressCallback.hpp"
#include "QAP_query.hpp"
#include "Rank1DSL.hpp"
#include "WindowExp.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Key pair: proving and verification
//

template <typename PAIRING>
class PPZK_Keypair
{
    typedef typename PAIRING::Fr Fr;
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;

public:
    PPZK_Keypair() = default;

    // only used for roundtrip marshalling tests
    PPZK_Keypair(const PPZK_ProvingKey<PAIRING>& pk,
                 const PPZK_VerificationKey<PAIRING> vk)
        : m_pk(pk),
          m_vk(vk)
    {}

    template <template <typename> class SYS>
    PPZK_Keypair(const SYS<Fr>& constraintSystem,
                 const std::size_t numCircuitInputs,
                 const PPZK_LagrangePoint<Fr>& lagrangeRand,
                 const PPZK_BlindGreeks<Fr, Fr>& blindRand,
                 ProgressCallback* callback = nullptr)
    {
        ProgressCallback_NOP<PAIRING> dummyNOP;
        ProgressCallback* dummy = callback ? callback : std::addressof(dummyNOP);
        dummy->majorSteps(8);

        // randomness
        const auto
            &point = lagrangeRand.point(),
            &alphaA = blindRand.alphaA(),
            &alphaB = blindRand.alphaB(),
            &alphaC = blindRand.alphaC(),
            &rA = blindRand.rA(),
            &rB = blindRand.rB(),
            &rC = blindRand.rC(),
            &beta = blindRand.beta(),
            &gamma = blindRand.gamma();

        const QAP_SystemPoint<SYS, Fr> qap(constraintSystem, numCircuitInputs, point);

        // ABCH
        QAP_QueryABC<SYS, Fr> ABCt(qap); // changed by QAP_QueryIC side-effect
        const QAP_QueryH<SYS, Fr> Ht(qap);

        // step 8 - G1 window table
        dummy->major(true);
        const WindowExp<G1> g1_table(g1_exp_count(qap, ABCt, Ht), callback);

        // step 7 - G2 window table
        dummy->major(true);
        const WindowExp<G2> g2_table(g2_exp_count(ABCt), callback);

        // step 6 - K
        dummy->major(true);
        auto K = ppzk_query_HK<PAIRING>(qap_query_K(qap, ABCt, rA, rB, beta), g1_table, callback);
#ifdef USE_ADD_SPECIAL
        batchSpecial(K);
#endif

        // step 5 - input consistency (side-effect: modifies ABCt query vector A)
        dummy->major(true);
        PPZK_QueryIC<PAIRING> ppzkIC(qap_query_IC(qap, ABCt, rA));
        ppzkIC.accumTable(g1_table, callback);

        // step 4 - A
        dummy->major(true);
        auto A = ppzk_query_ABC(ABCt.vecA(), rA, alphaA, g1_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        batchSpecial(A);
#endif

        // step 3 - B
        dummy->major(true);
        auto B = ppzk_query_ABC(ABCt.vecB(), rB, alphaB, g2_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        batchSpecial(B);
#endif

        // step 2 - C
        dummy->major(true);
        auto C = ppzk_query_ABC(ABCt.vecC(), rC, alphaC, g1_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        batchSpecial(C);
#endif

        // step 1 - H
        dummy->major(true);
        auto H = ppzk_query_HK<PAIRING>(Ht.vec(), g1_table, callback);

        m_pk = PPZK_ProvingKey<PAIRING>(std::move(A),
                                        std::move(B),
                                        std::move(C),
                                        std::move(H),
                                        std::move(K));

        m_vk = PPZK_VerificationKey<PAIRING>(alphaA * G2::one(),
                                             alphaB * G1::one(),
                                             alphaC * G2::one(),
                                             gamma * G2::one(),
                                             (gamma * beta) * G1::one(),
                                             (gamma * beta) * G2::one(),
                                             (rC * qap.compute_Z()) * G2::one(),
                                             std::move(ppzkIC));
    }

    const PPZK_ProvingKey<PAIRING>& pk() const { return m_pk; }
    const PPZK_VerificationKey<PAIRING>& vk() const { return m_vk; }

    bool operator== (const PPZK_Keypair& other) const {
        return
            pk() == other.pk() &&
            vk() == other.vk();
    }

    bool operator!= (const PPZK_Keypair& other) const {
        return ! (*this == other);
    }

    void marshal_out(std::ostream& os) const {
        pk().marshal_out(os);
        vk().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_pk.marshal_in(is) &&
            m_vk.marshal_in(is);
    }

    void clear() {
        m_pk.clear();
        m_vk.clear();
    }

    bool empty() const {
        return
            m_pk.empty() ||
            m_vk.empty();
    }

private:
    PPZK_ProvingKey<PAIRING> m_pk;
    PPZK_VerificationKey<PAIRING> m_vk;
};

template <typename PAIRING>
std::ostream& operator<< (std::ostream& os, const PPZK_Keypair<PAIRING>& a) {
    a.marshal_out(os);
    return os;
}

template <typename PAIRING>
std::istream& operator>> (std::istream& is, PPZK_Keypair<PAIRING>& a) {
    if (! a.marshal_in(is)) a.clear();
    return is;
}

} // namespace snarklib

#endif
