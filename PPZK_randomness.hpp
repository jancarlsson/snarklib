#ifndef _SNARKLIB_PPZK_RANDOMNESS_HPP_
#define _SNARKLIB_PPZK_RANDOMNESS_HPP_

#include <cassert>
#include <cstdint>
#include <istream>
#include <iostream>
#include <ostream>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// proving and verification keypair randomness
// (SECURITY CAUTION: destroy unblinded entropy after keypair generation)
//

// If BLIND = FR, then randomness is exposed in the clear (dangerous).
// If BLIND = Pairing<PAIRING::G1, PAIRING::G2>, then randomness is blinded (safe).
template <typename FR, typename BLIND>
class PPZK_KeypairRandomness
{
public:
    PPZK_KeypairRandomness() = default;

    PPZK_KeypairRandomness(const std::size_t degree)
        : m_point(degree + 1, BLIND::zero()),
          m_alphaA(FR::random() * BLIND::one()),
          m_alphaB(FR::random() * BLIND::one()),
          m_alphaC(FR::random() * BLIND::one()),
          m_beta(FR::random() * BLIND::one()),
          m_gamma(FR::random() * BLIND::one())
    {
#ifdef USE_ASSERT
        assert(degree > 0);
#endif

        // point powers from 0 to degree inclusive
        const auto point = FR::random();
        auto t = FR::one();
        for (auto& u : m_point) {
            u = t * BLIND::one();
            t *= point;
        }

        // rC = rA * rB
        const auto a = FR::random();
        m_rA = a * BLIND::one();
        m_rB = FR::random() * BLIND::one();
        m_rC = a * m_rB;
    }

    std::size_t degree() const { return m_point.size() - 1; }
    const std::vector<BLIND>& pointvec() const { return m_point; }
    const BLIND& point() const { return pointvec()[1]; }

    const BLIND& alphaA() const { return m_alphaA; }
    const BLIND& alphaB() const { return m_alphaB; }
    const BLIND& alphaC() const { return m_alphaC; }

    const BLIND& rA() const { return m_rA; }
    const BLIND& rB() const { return m_rB; }
    const BLIND& rC() const { return m_rC; }

    const BLIND& beta() const { return m_beta; }
    const BLIND& gamma() const { return m_gamma; }

    // SECURITY CAUTION: destination of stream is trusted
    void marshal_out(std::ostream& os) const {
        snarklib::marshal_out(os, pointvec());
        alphaA().marshal_out(os);
        alphaB().marshal_out(os);
        alphaC().marshal_out(os);
        rA().marshal_out(os);
        rB().marshal_out(os);
        rC().marshal_out(os);
        beta().marshal_out(os);
        gamma().marshal_out(os);
    }

    // SECURITY CAUTION: source of stream is trusted
    bool marshal_in(std::istream& is) {
        return
            snarklib::marshal_in(is, m_point) &&
            m_alphaA.marshal_in(is) &&
            m_alphaB.marshal_in(is) &&
            m_alphaC.marshal_in(is) &&
            m_rA.marshal_in(is) &&
            m_rB.marshal_in(is) &&
            m_rC.marshal_in(is) &&
            m_beta.marshal_in(is) &&
            m_gamma.marshal_in(is);
    }

    void clear() {
        m_point.clear();
        m_alphaA.clear();
        m_alphaB.clear();
        m_alphaC.clear();
        m_rA.clear();
        m_rB.clear();
        m_rC.clear();
        m_beta.clear();
        m_gamma.clear();
    }

    bool empty() const {
        return
            m_point.empty() ||
            m_alphaA.isZero() ||
            m_alphaB.isZero() ||
            m_alphaC.isZero() ||
            m_rA.isZero() ||
            m_rB.isZero() ||
            m_rC.isZero() ||
            m_beta.isZero() ||
            m_gamma.isZero();
    }

private:
    std::vector<BLIND> m_point;
    BLIND m_alphaA, m_alphaB, m_alphaC;
    BLIND m_rA, m_rB, m_rC;
    BLIND m_beta, m_gamma;
};

template <typename FR, typename BLIND>
std::ostream& operator<< (std::ostream& os,
                          const PPZK_KeypairRandomness<FR, BLIND>& a) {
    a.marshal_out(os);
    return os;
}

template <typename FR, typename BLIND>
std::istream& operator>> (std::istream& is,
                          PPZK_KeypairRandomness<FR, BLIND>& a) {
    if (! a.marshal_in(is)) a.clear();
    return is;
}

////////////////////////////////////////////////////////////////////////////////
// proof generation randomness
//

template <typename FR>
class PPZK_ProofRandomness
{
public:
    PPZK_ProofRandomness() = default;

    PPZK_ProofRandomness(const int dummy)
        : m_d1(FR::random()),
          m_d2(FR::random()),
          m_d3(FR::random())
    {}

    const FR& d1() const { return m_d1; }
    const FR& d2() const { return m_d2; }
    const FR& d3() const { return m_d3; }

    void marshal_out(std::ostream& os) const {
        d1().marshal_out(os);
        d2().marshal_out(os);
        d3().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_d1.marshal_in(is) &&
            m_d2.marshal_in(is) &&
            m_d3.marshal_in(is);
    }

    void clear() {
        m_d1.clear();
        m_d2.clear();
        m_d3.clear();
    }

    bool empty() const {
        return
            m_d1.isZero() ||
            m_d2.isZero() ||
            m_d3.isZero();
    }

private:
    FR m_d1, m_d2, m_d3;
};

template <typename FR>
std::ostream& operator<< (std::ostream& os, const PPZK_ProofRandomness<FR>& a) {
    a.marshal_out(os);
    return os;
}

template <typename FR>
std::istream& operator>> (std::istream& is, PPZK_ProofRandomness<FR>& a) {
    if (! a.marshal_in(is)) a.clear();
    return is;
}

} // namespace snarklib

#endif
