#ifndef _SNARKLIB_EC_MNT6_INIT_GROUPS_HPP_
#define _SNARKLIB_EC_MNT6_INIT_GROUPS_HPP_

#include <snarklib/EC.hpp>
#include <snarklib/EC_MNT6_GroupCurve.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT6
// initialize groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class MNT6_InitGroups : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef MNT6_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

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
        G1::params.G_zero(Fq::zero(), Fq::one(), Fq::zero());
        G1::params.G_one(
            "336685752883082228109289846353937104185698209371404178342968838739115829740084426881123453",
            "402596290139780989709332707716568920777622032073762749862342374583908837063963736098549800",
            Fq::one());

        G1::params.wnaf_window_table_clear();
        G1::params.wnaf_window_table(11);
        G1::params.wnaf_window_table(24);
        G1::params.wnaf_window_table(60);
        G1::params.wnaf_window_table(127);

        G1::params.fixed_base_exp_window_table_clear();
        G1::params.fixed_base_exp_window_table(1);
        G1::params.fixed_base_exp_window_table(4);
        G1::params.fixed_base_exp_window_table(10);
        G1::params.fixed_base_exp_window_table(25);
        G1::params.fixed_base_exp_window_table(60);
        G1::params.fixed_base_exp_window_table(146);
        G1::params.fixed_base_exp_window_table(350);
        G1::params.fixed_base_exp_window_table(845);
        G1::params.fixed_base_exp_window_table(1840);
        G1::params.fixed_base_exp_window_table(3904);
        G1::params.fixed_base_exp_window_table(11309);
        G1::params.fixed_base_exp_window_table(24016);
        G1::params.fixed_base_exp_window_table(72289);
        G1::params.fixed_base_exp_window_table(138413);
        G1::params.fixed_base_exp_window_table(156390);
        G1::params.fixed_base_exp_window_table(562560);
        G1::params.fixed_base_exp_window_table(1036742);
        G1::params.fixed_base_exp_window_table(2053819);
        G1::params.fixed_base_exp_window_table(4370224);
        G1::params.fixed_base_exp_window_table(8215704);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(42682375);

        G2::params.G_zero(Fq3::zero(), Fq3::one(), Fq3::zero());
        G2::params.G_one(
            Fq3("421456435772811846256826561593908322288509115489119907560382401870203318738334702321297427",
                "103072927438548502463527009961344915021167584706439945404959058962657261178393635706405114",
                "143029172143731852627002926324735183809768363301149009204849580478324784395590388826052558"),
            Fq3("464673596668689463130099227575639512541218133445388869383893594087634649237515554342751377",
                "100642907501977375184575075967118071807821117960152743335603284583254620685343989304941678",
                "123019855502969896026940545715841181300275180157288044663051565390506010149881373807142903"),
            Fq3::one());

        G2::params.wnaf_window_table_clear();
        G2::params.wnaf_window_table(5);
        G2::params.wnaf_window_table(15);
        G2::params.wnaf_window_table(39);
        G2::params.wnaf_window_table(109);

        G2::params.fixed_base_exp_window_table_clear();
        G2::params.fixed_base_exp_window_table(1);
        G2::params.fixed_base_exp_window_table(4);
        G2::params.fixed_base_exp_window_table(10);
        G2::params.fixed_base_exp_window_table(25);
        G2::params.fixed_base_exp_window_table(60);
        G2::params.fixed_base_exp_window_table(144);
        G2::params.fixed_base_exp_window_table(346);
        G2::params.fixed_base_exp_window_table(819);
        G2::params.fixed_base_exp_window_table(1782);
        G2::params.fixed_base_exp_window_table(4002);
        G2::params.fixed_base_exp_window_table(10870);
        G2::params.fixed_base_exp_window_table(18023);
        G2::params.fixed_base_exp_window_table(43161);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(149743);
        G2::params.fixed_base_exp_window_table(551844);
        G2::params.fixed_base_exp_window_table(1041828);
        G2::params.fixed_base_exp_window_table(1977372);
        G2::params.fixed_base_exp_window_table(3703620);
        G2::params.fixed_base_exp_window_table(7057237);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(38554492);
    }
};

} // namespace snarklib

#endif
