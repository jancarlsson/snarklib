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
    static PPZK_KeypairRandomness<Fr> randomness() {
        return PPZK_KeypairRandomness<Fr>(0);
    }

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
                 const PPZK_KeypairRandomness<Fr>& keyRand,
                 ProgressCallback* callback = nullptr)
    {
        ProgressCallback_NOP<PAIRING> dummyNOP;
        ProgressCallback* dummy = callback ? callback : std::addressof(dummyNOP);
        dummy->majorSteps(8);

        // randomness
        const auto
            &point = keyRand.point(),
            &alphaA = keyRand.alphaA(),
            &alphaB = keyRand.alphaB(),
            &alphaC = keyRand.alphaC(),
            &rA = keyRand.rA(),
            &rB = keyRand.rB(),
            &rC = keyRand.rC(),
            &beta = keyRand.beta(),
            &gamma = keyRand.gamma();

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
        const QAP_QueryK<SYS, Fr> Kt(qap, ABCt, rA, rB, beta);
        const BlockVector<Fr> Ktb(BlockVector<Fr>::space(Kt.vec()), 0, Kt.vec());
        PPZK_QueryK<PAIRING> Kp(Ktb);
        Kp.accumTable(g1_table, callback);
#ifdef USE_ADD_SPECIAL
        Kp.batchSpecial();
#endif

        // side-effect: this modifies ABCt query vector A
        const QAP_QueryIC<SYS, Fr> qapIC(qap, ABCt, rA);

        // step 5 - A
        dummy->major(true);
        const BlockVector<Fr> Atb(BlockVector<Fr>::space(ABCt.vecA()), 0, ABCt.vecA());
        PPZK_QueryA<PAIRING> Ap(Atb, rA, alphaA);
        Ap.accumTable(g1_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        Ap.batchSpecial();
#endif

        // step 4 - B
        dummy->major(true);
        const BlockVector<Fr> Btb(BlockVector<Fr>::space(ABCt.vecB()), 0, ABCt.vecB());
        PPZK_QueryB<PAIRING> Bp(Btb, rB, alphaB);
        Bp.accumTable(g2_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        Bp.batchSpecial();
#endif

        // step 3 - C
        dummy->major(true);
        const BlockVector<Fr> Ctb(BlockVector<Fr>::space(ABCt.vecC()), 0, ABCt.vecC());
        PPZK_QueryC<PAIRING> Cp(Ctb, rC, alphaC);
        Cp.accumTable(g1_table, g1_table, callback);
#ifdef USE_ADD_SPECIAL
        Cp.batchSpecial();
#endif

        // step 2 - H
        dummy->major(true);
        const BlockVector<Fr> Htb(BlockVector<Fr>::space(Ht.vec()), 0, Ht.vec());
        PPZK_QueryH<PAIRING> Hp(Htb);
        Hp.accumTable(g1_table, callback);

        m_pk = PPZK_ProvingKey<PAIRING>(Ap.vec(),
                                        Bp.vec(),
                                        Cp.vec(),
                                        Hp.vvec(),
                                        Kp.vvec());

        // step 1 - input consistency
        dummy->major(true);
        PPZK_QueryIC<PAIRING> ppzkIC(qapIC.vec());
        ppzkIC.accumTable(g1_table, callback);

        m_vk = PPZK_VerificationKey<PAIRING>(alphaA * G2::one(),
                                             alphaB * G1::one(),
                                             alphaC * G2::one(),
                                             gamma * G2::one(),
                                             (gamma * beta) * G1::one(),
                                             (gamma * beta) * G2::one(),
                                             (rC * qap.compute_Z()) * G2::one(),
                                             ppzkIC);
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
