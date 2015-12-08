#ifndef _SNARKLIB_EC_MNT4_INIT_GROUPS_HPP_
#define _SNARKLIB_EC_MNT4_INIT_GROUPS_HPP_

#include <snarklib/EC.hpp>
#include <snarklib/EC_MNT4_GroupCurve.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT4
// initialize groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class MNT4_InitGroups : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef MNT4_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

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
        G1::params.G_one(
            "60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838",
            "363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306",
            Fq::one());

        G1::params.wnaf_window_table_clear();
        G1::params.wnaf_window_table(11);
        G1::params.wnaf_window_table(24);
        G1::params.wnaf_window_table(60);
        G1::params.wnaf_window_table(127);

        G1::params.fixed_base_exp_window_table_clear();
        G1::params.fixed_base_exp_window_table(1);
        G1::params.fixed_base_exp_window_table(5);
        G1::params.fixed_base_exp_window_table(10);
        G1::params.fixed_base_exp_window_table(25);
        G1::params.fixed_base_exp_window_table(60);
        G1::params.fixed_base_exp_window_table(144);
        G1::params.fixed_base_exp_window_table(345);
        G1::params.fixed_base_exp_window_table(855);
        G1::params.fixed_base_exp_window_table(1805);
        G1::params.fixed_base_exp_window_table(3912);
        G1::params.fixed_base_exp_window_table(11265);
        G1::params.fixed_base_exp_window_table(27898);
        G1::params.fixed_base_exp_window_table(57597);
        G1::params.fixed_base_exp_window_table(145299);
        G1::params.fixed_base_exp_window_table(157205);
        G1::params.fixed_base_exp_window_table(601601);
        G1::params.fixed_base_exp_window_table(1107377);
        G1::params.fixed_base_exp_window_table(1789647);
        G1::params.fixed_base_exp_window_table(4392627);
        G1::params.fixed_base_exp_window_table(8221211);
        G1::params.fixed_base_exp_window_table(0);
        G1::params.fixed_base_exp_window_table(42363731);

        G2::params.G_zero(Fq2::zero(), Fq2::one(), Fq2::zero());
        G2::params.G_one(
            Fq2("438374926219350099854919100077809681842783509163790991847867546339851681564223481322252708",
                "37620953615500480110935514360923278605464476459712393277679280819942849043649216370485641"),
            Fq2("37437409008528968268352521034936931842973546441370663118543015118291998305624025037512482",
                "424621479598893882672393190337420680597584695892317197646113820787463109735345923009077489"),
            Fq2::one());

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
        G2::params.fixed_base_exp_window_table(143);
        G2::params.fixed_base_exp_window_table(345);
        G2::params.fixed_base_exp_window_table(821);
        G2::params.fixed_base_exp_window_table(1794);
        G2::params.fixed_base_exp_window_table(3920);
        G2::params.fixed_base_exp_window_table(11301);
        G2::params.fixed_base_exp_window_table(18960);
        G2::params.fixed_base_exp_window_table(44199);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(150800);
        G2::params.fixed_base_exp_window_table(548695);
        G2::params.fixed_base_exp_window_table(1051769);
        G2::params.fixed_base_exp_window_table(2023926);
        G2::params.fixed_base_exp_window_table(3787109);
        G2::params.fixed_base_exp_window_table(7107480);
        G2::params.fixed_base_exp_window_table(0);
        G2::params.fixed_base_exp_window_table(38760027);
    }
};

} // namespace snarklib

#endif
