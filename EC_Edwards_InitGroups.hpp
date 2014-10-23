#ifndef _SNARKLIB_EC_EDWARDS_INIT_GROUPS_HPP_
#define _SNARKLIB_EC_EDWARDS_INIT_GROUPS_HPP_

#include "EC.hpp"
#include "EC_Edwards_GroupCurve.hpp"
#include "Group.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Edwards (80 bits)
// initialize groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class Edwards_InitGroups : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef Edwards_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

public:
    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq3 Fq3; // twist field for G2

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq3, Fr, CURVE> G2;

    // initialize group parameters
    static void initParams()
    {
        G1::params.G_zero(Fq::zero(), Fq::one());
        G1::params.G_one("3713709671941291996998665608188072510389821008693530490",
                         "4869953702976555123067178261685365085639705297852816679");

        G1::params.wnaf_window_table_clear();
        G1::params.wnaf_window_table(9);
        G1::params.wnaf_window_table(14);
        G1::params.wnaf_window_table(24);
        G1::params.wnaf_window_table(117);

        G1::params.fixed_base_exp_window_table_clear();
        G1::params.fixed_base_exp_window_table(1);
        G1::params.fixed_base_exp_window_table(4);
        G1::params.fixed_base_exp_window_table(10);
        G1::params.fixed_base_exp_window_table(25);
        G1::params.fixed_base_exp_window_table(60);
        G1::params.fixed_base_exp_window_table(149);
        G1::params.fixed_base_exp_window_table(370);
        G1::params.fixed_base_exp_window_table(849);
        G1::params.fixed_base_exp_window_table(1765);
        G1::params.fixed_base_exp_window_table(4430);
        G1::params.fixed_base_exp_window_table(13389);
        G1::params.fixed_base_exp_window_table(15368);
        G1::params.fixed_base_exp_window_table(74912);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(438107);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(1045626);
        G1::params.fixed_base_exp_window_table(1577434);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(17350594);
        G1::params.fixed_base_exp_window_table(0);

        G2::params.G_zero(Fq3::zero(), Fq3::one());
        G2::params.G_one(Fq3("4531683359223370252210990718516622098304721701253228128",
                             "5339624155305731263217400504407647531329993548123477368",
                             "3964037981777308726208525982198654699800283729988686552"),
                         Fq3("364634864866983740775341816274081071386963546650700569",
                             "3264380230116139014996291397901297105159834497864380415",
                             "3504781284999684163274269077749440837914479176282903747"));

        G2::params.wnaf_window_table_clear();
        G2::params.wnaf_window_table(6);
        G2::params.wnaf_window_table(12);
        G2::params.wnaf_window_table(42);
        G2::params.wnaf_window_table(97);

        G2::params.fixed_base_exp_window_table_clear();
        G2::params.fixed_base_exp_window_table(1);
        G2::params.fixed_base_exp_window_table(5);
        G2::params.fixed_base_exp_window_table(11);
        G2::params.fixed_base_exp_window_table(26);
        G2::params.fixed_base_exp_window_table(61);
        G2::params.fixed_base_exp_window_table(146);
        G2::params.fixed_base_exp_window_table(357);
        G2::params.fixed_base_exp_window_table(823);
        G2::params.fixed_base_exp_window_table(1589);
        G2::params.fixed_base_exp_window_table(4136);
        G2::params.fixed_base_exp_window_table(14298);
        G2::params.fixed_base_exp_window_table(16745);
        G2::params.fixed_base_exp_window_table(51769);
        G2::params.fixed_base_exp_window_table(99811);
        G2::params.fixed_base_exp_window_table(193307);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(907185);
        G2::params.fixed_base_exp_window_table(1389683);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(6752696);
        G2::params.fixed_base_exp_window_table(193642895);
        G2::params.fixed_base_exp_window_table(226760202);
    }
};

} // namespace snarklib

#endif
