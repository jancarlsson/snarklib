#ifndef _SNARKLIB_PPZK_RANDOMNESS_HPP_
#define _SNARKLIB_PPZK_RANDOMNESS_HPP_

#include <istream>
#include <ostream>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// proving and verification keypair randomness
// (SECURITY CAUTION: destroy after keypair generation)
//

template <typename FR>
class PPZK_KeypairRandomness
{
public:
    PPZK_KeypairRandomness() = default;

    PPZK_KeypairRandomness(const int dummy)
        : m_point(FR::random()),
          m_alphaA(FR::random()),
          m_alphaB(FR::random()),
          m_alphaC(FR::random()),
          m_rA(FR::random()),
          m_rB(FR::random()),
          m_rC(m_rA * m_rB),
          m_beta(FR::random()),
          m_gamma(FR::random())
    {}

    const FR& point() const { return m_point; }

    const FR& alphaA() const { return m_alphaA; }
    const FR& alphaB() const { return m_alphaB; }
    const FR& alphaC() const { return m_alphaC; }

    const FR& rA() const { return m_rA; }
    const FR& rB() const { return m_rB; }
    const FR& rC() const { return m_rC; }

    const FR& beta() const { return m_beta; }
    const FR& gamma() const { return m_gamma; }

    // SECURITY CAUTION: destination of stream is trusted
    void marshal_out(std::ostream& os) const {
        point().marshal_out(os);
        alphaA().marshal_out(os);
        alphaB().marshal_out(os);
        alphaC().marshal_out(os);
        rA().marshal_out(os);
        rB().marshal_out(os);
        beta().marshal_out(os);
        gamma().marshal_out(os);
    }

    // SECURITY CAUTION: source of stream is trusted
    bool marshal_in(std::istream& is) {
        const bool rc =
            m_point.marshal_in(is) &&
            m_alphaA.marshal_in(is) &&
            m_alphaB.marshal_in(is) &&
            m_alphaC.marshal_in(is) &&
            m_rA.marshal_in(is) &&
            m_rB.marshal_in(is) &&
            m_beta.marshal_in(is) &&
            m_gamma.marshal_in(is);

        m_rC = m_rA * m_rB;

        return rc;
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
            m_point.isZero() ||
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
    FR m_point;
    FR m_alphaA, m_alphaB, m_alphaC;
    FR m_rA;
    FR m_rB;
    FR m_rC;
    FR m_beta, m_gamma;
};

template <typename FR>
std::ostream& operator<< (std::ostream& os, const PPZK_KeypairRandomness<FR>& a) {
    a.marshal_out(os);
    return os;
}

template <typename FR>
std::istream& operator>> (std::istream& is, PPZK_KeypairRandomness<FR>& a) {
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
