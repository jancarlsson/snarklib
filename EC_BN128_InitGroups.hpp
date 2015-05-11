#ifndef _SNARKLIB_EC_BN128_INIT_GROUPS_HPP_
#define _SNARKLIB_EC_BN128_INIT_GROUPS_HPP_

#include <snarklib/EC.hpp>
#include <snarklib/EC_BN128_GroupCurve.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Barreto-Naehrig (128 bits)
// initialize groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class BN128_InitGroups : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef BN128_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

public:
    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq2 Fq2; // twist field for G2

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq2, Fr, CURVE> G2;

    // initialize group parameters
    static void initParams()
    {
        G1::params.G_zero(Fq::zero(), Fq::one(), Fq::zero());
        G1::params.G_one("1", "2", Fq::one());

        G1::params.wnaf_window_table_clear();
        G1::params.wnaf_window_table(11);
        G1::params.wnaf_window_table(24);
        G1::params.wnaf_window_table(60);
        G1::params.wnaf_window_table(127);

        G1::params.fixed_base_exp_window_table_clear();
        G1::params.fixed_base_exp_window_table(1);
        G1::params.fixed_base_exp_window_table(5);
        G1::params.fixed_base_exp_window_table(11);
        G1::params.fixed_base_exp_window_table(32);
        G1::params.fixed_base_exp_window_table(55);
        G1::params.fixed_base_exp_window_table(162);
        G1::params.fixed_base_exp_window_table(360);
        G1::params.fixed_base_exp_window_table(815);
        G1::params.fixed_base_exp_window_table(2373);
        G1::params.fixed_base_exp_window_table(6978);
        G1::params.fixed_base_exp_window_table(7122);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(57818);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(169679);
        G1::params.fixed_base_exp_window_table(439759);
        G1::params.fixed_base_exp_window_table(936073);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(4666555);
        G1::params.fixed_base_exp_window_table(7580404);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(34552892);

        G2::params.G_zero(Fq2::zero(), Fq2::one(), Fq2::zero());
        G2::params.G_one(
            Fq2("10857046999023057135944570762232829481370756359578518086990519993285655852781",
                "11559732032986387107991004021392285783925812861821192530917403151452391805634"),
            Fq2("8495653923123431417604973247489272438418190587263600148770280649306958101930",
                "4082367875863433681332203403145435568316851327593401208105741076214120093531"),
            Fq2::one());

        G2::params.wnaf_window_table_clear();
        G2::params.wnaf_window_table(5);
        G2::params.wnaf_window_table(15);
        G2::params.wnaf_window_table(39);
        G2::params.wnaf_window_table(109);

        G2::params.fixed_base_exp_window_table_clear();
        G2::params.fixed_base_exp_window_table(1);
        G2::params.fixed_base_exp_window_table(5);
        G2::params.fixed_base_exp_window_table(10);
        G2::params.fixed_base_exp_window_table(25);
        G2::params.fixed_base_exp_window_table(59);
        G2::params.fixed_base_exp_window_table(154);
        G2::params.fixed_base_exp_window_table(334);
        G2::params.fixed_base_exp_window_table(743);
        G2::params.fixed_base_exp_window_table(2034);
        G2::params.fixed_base_exp_window_table(4988);
        G2::params.fixed_base_exp_window_table(8888);
        G2::params.fixed_base_exp_window_table(26271);
        G2::params.fixed_base_exp_window_table(39768);
        G2::params.fixed_base_exp_window_table(106276);
        G2::params.fixed_base_exp_window_table(141703);
        G2::params.fixed_base_exp_window_table(462423);
        G2::params.fixed_base_exp_window_table(926872);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(4873049);
        G2::params.fixed_base_exp_window_table(5706708);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(31673815);
    }
};

} // namespace snarklib

#endif
