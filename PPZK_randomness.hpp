#ifndef _SNARKLIB_PPZK_RANDOMNESS_HPP_
#define _SNARKLIB_PPZK_RANDOMNESS_HPP_

#include <cstdint>
#include <istream>
#include <iostream>
#include <ostream>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// proving and verification keypair randomness
// (SECURITY CAUTION: destroy unblinded entropy after keypair generation)
//

// not blinded
template <typename FR>
class PPZK_LagrangePoint
{
public:
    PPZK_LagrangePoint() = default;

    PPZK_LagrangePoint(const int dummy)
        : m_point(FR::random())
    {}

    const FR& point() const {
        return m_point;
    }

    void marshal_out(std::ostream& os) const {
        point().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return m_point.marshal_in(is);
    }

    void clear() {
        m_point.clear();
    }

    bool empty() const {
        return point().isZero();
    }

private:
    FR m_point;
};

template <typename FR>
std::ostream& operator<< (std::ostream& os,
                          const PPZK_LagrangePoint<FR>& a) {
    a.marshal_out(os);
    return os;
}

template <typename FR>
std::istream& operator>> (std::istream& is,
                          PPZK_LagrangePoint<FR>& a) {
    if (! a.marshal_in(is)) a.clear();
    return is;
}

// If BLIND = FR, then randomness is exposed in the clear (dangerous).
// If BLIND = Pairing<PAIRING::G1, PAIRING::G2>, then randomness is blinded (safe).
template <typename FR, typename BLIND>
class PPZK_BlindGreeks
{
public:
    PPZK_BlindGreeks() = default;

    PPZK_BlindGreeks(const int dummy)
        : m_rB(FR::random() * BLIND::one()),
          m_gamma(FR::random() * BLIND::one())
    {
        const auto
            alphaA = FR::random(),
            alphaB = FR::random(),
            alphaC = FR::random(),
            rA = FR::random(),
            beta = FR::random();

        m_alphaA = alphaA * BLIND::one();
        m_alphaB = alphaB * BLIND::one();
        m_alphaC = alphaC * BLIND::one();

        m_rA = rA * BLIND::one();
        m_rC = rA * m_rB;

        m_beta = beta * BLIND::one();

        m_alphaA_rA = alphaA * m_rA;
        m_alphaB_rB = alphaB * m_rB;
        m_alphaC_rC = alphaC * m_rC;

        m_beta_rA = beta * m_rA;
        m_beta_rB = beta * m_rB;
        m_beta_rC = beta * m_rC;

        m_beta_gamma = beta * m_gamma;
    }

    const BLIND& alphaA() const { return m_alphaA; }
    const BLIND& alphaB() const { return m_alphaB; }
    const BLIND& alphaC() const { return m_alphaC; }

    const BLIND& rA() const { return m_rA; }
    const BLIND& rB() const { return m_rB; }
    const BLIND& rC() const { return m_rC; }

    const BLIND& beta() const { return m_beta; }
    const BLIND& gamma() const { return m_gamma; }

    const BLIND& alphaA_rA() const { return m_alphaA_rA; }
    const BLIND& alphaB_rB() const { return m_alphaB_rB; }
    const BLIND& alphaC_rC() const { return m_alphaC_rC; }

    const BLIND& beta_rA() const { return m_beta_rA; }
    const BLIND& beta_rB() const { return m_beta_rB; }
    const BLIND& beta_rC() const { return m_beta_rC; }

    const BLIND& beta_gamma() const { return m_beta_gamma; }

    void marshal_out(std::ostream& os) const {
        alphaA().marshal_out(os);
        alphaB().marshal_out(os);
        alphaC().marshal_out(os);

        rA().marshal_out(os);
        rB().marshal_out(os);
        rC().marshal_out(os);

        beta().marshal_out(os);
        gamma().marshal_out(os);

        alphaA_rA().marshal_out(os);
        alphaB_rB().marshal_out(os);
        alphaC_rC().marshal_out(os);

        beta_rA().marshal_out(os);
        beta_rB().marshal_out(os);
        beta_rC().marshal_out(os);

        beta_gamma().marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        return
            m_alphaA.marshal_in(is) &&
            m_alphaB.marshal_in(is) &&
            m_alphaC.marshal_in(is) &&

            m_rA.marshal_in(is) &&
            m_rB.marshal_in(is) &&
            m_rC.marshal_in(is) &&

            m_beta.marshal_in(is) &&
            m_gamma.marshal_in(is) &&

            m_alphaA_rA.marshal_in(is) &&
            m_alphaB_rB.marshal_in(is) &&
            m_alphaC_rC.marshal_in(is) &&

            m_beta_rA.marshal_in(is) &&
            m_beta_rB.marshal_in(is) &&
            m_beta_rC.marshal_in(is) &&

            m_beta_gamma.marshal_in(is);
    }

    void clear() {
        m_alphaA.clear();
        m_alphaB.clear();
        m_alphaC.clear();

        m_rA.clear();
        m_rB.clear();
        m_rC.clear();

        m_beta.clear();
        m_gamma.clear();

        m_alphaA_rA.clear();
        m_alphaB_rB.clear();
        m_alphaC_rC.clear();

        m_beta_rA.clear();
        m_beta_rB.clear();
        m_beta_rC.clear();

        m_beta_gamma.clear();
    }

    bool empty() const {
        return
            alphaA().isZero() ||
            alphaB().isZero() ||
            alphaC().isZero() ||

            rA().isZero() ||
            rB().isZero() ||
            rC().isZero() ||

            beta().isZero() ||
            gamma().isZero() ||

            alphaA_rA().isZero() ||
            alphaB_rB().isZero() ||
            alphaC_rC().isZero() ||

            beta_rA().isZero() ||
            beta_rB().isZero() ||
            beta_rC().isZero() ||

            beta_gamma().isZero();
    }

private:
    // random parameters
    BLIND m_alphaA, m_alphaB, m_alphaC;
    BLIND m_rA, m_rB, m_rC;
    BLIND m_beta, m_gamma;

    // products of random parameters
    BLIND m_alphaA_rA, m_alphaB_rB, m_alphaC_rC;
    BLIND m_beta_rA, m_beta_rB, m_beta_rC;
    BLIND m_beta_gamma;
};

template <typename FR, typename BLIND>
std::ostream& operator<< (std::ostream& os,
                          const PPZK_BlindGreeks<FR, BLIND>& a) {
    a.marshal_out(os);
    return os;
}

template <typename FR, typename BLIND>
std::istream& operator>> (std::istream& is,
                          PPZK_BlindGreeks<FR, BLIND>& a) {
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
