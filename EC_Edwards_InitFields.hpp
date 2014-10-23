#ifndef _SNARKLIB_EC_EDWARDS_INIT_FIELDS_HPP_
#define _SNARKLIB_EC_EDWARDS_INIT_FIELDS_HPP_

#include "EC.hpp"
#include "EC_Edwards_Modulus.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Edwards (80 bits)
// initialize fields
//

// MODULUS is: modulus_r; modulus_q
template <mp_size_t N, const BigInt<N>& MODULUS>
class Edwards_InitFields : public ECInitField<N, MODULUS, Edwards_InitFields<N, MODULUS>> // CRTP
{
    typedef ECInitField<N, MODULUS, Edwards_InitFields<N, MODULUS>> BASE;

public:
    static const BigInt<Edwards_Modulus::r_limbs>& modulus_r() {
        return Edwards_Modulus::modulus_r();
    }

    static const BigInt<Edwards_Modulus::q_limbs>& modulus_q() {
        return Edwards_Modulus::modulus_q();
    }

    typedef typename BASE::Fp F;
    typedef typename BASE::Fp2 F2; // required by F6 cyclotomic_squared()
    typedef typename BASE::Fp3 F3;
    typedef typename BASE::Fp32 F6;

    static void initModulusR()
    {
        if (8 == sizeof(mp_limb_t)) {
            F::params.Rsquared("621738487827897760168419760282818735947979812540885779");
            F::params.Rcubed("899968968216802386013510389846941393831065658679774050");
            F::params.inv(0xdde553277fffffff);
        }

        if (4 == sizeof(mp_limb_t))
        {
            F::params.Rsquared("621738487827897760168419760282818735947979812540885779");
            F::params.Rcubed("899968968216802386013510389846941393831065658679774050");
            F::params.inv(0x7fffffff);
        }

        F::params.num_bits(181);
        F::params.s(31);
        F::params.t_minus_1_over_2("361472142418481002384052044093753675292693287");
        F::params.multiplicative_generator("19");
        F::params.root_of_unity("695314865466598274460565335217615316274564719601897184");
        F::params.nqr_to_t("1326707053668679463752768729767248251415639579872144553");
    }

    static void initModulusQ()
    {
        if (8 == sizeof(mp_limb_t)) {
            F::params.Rsquared("5943559676554581037560514598978484097352477055348195432");
            F::params.Rcubed("1081560488703514202058739223469726982199727506489234349");
            F::params.inv(0x76eb690b7fffffff);
        }

        if (4 == sizeof(mp_limb_t)) {
            F::params.Rsquared("5943559676554581037560514598978484097352477055348195432");
            F::params.Rcubed("1081560488703514202058739223469726982199727506489234349");
            F::params.inv(0x7fffffff);
        }

        F::params.num_bits(183);
        F::params.s(31);
        F::params.t_minus_1_over_2("1445888569673924009536208175329020776442194187");
        F::params.multiplicative_generator("61");
        F::params.root_of_unity("4692813029219384139894873043933463717810008194158530536");
        F::params.nqr_to_t("2626736066325740702418554487368721595489070118548299138");

        F3::params.s(31);
        F3::params.t_minus_1_over_2("55760183704072378092907654676152317178531104407263430256321995781805829544575551831417485592515824843119665712310518678891618803500033228219447095278582562");
        F3::params.non_residue("61");
        F3::params.nqr_to_t("104810943629412208121981114244673004633270996333237516", "0", "0");
        F3::params.Frobenius_coeffs_c1(0, "1");
        F3::params.Frobenius_coeffs_c1(1, "1073752683758513276629212192812154536507607213288832061");
        F3::params.Frobenius_coeffs_c1(2, "5136291436651207728317994048073823738016144056504959939");
        F3::params.Frobenius_coeffs_c2(0, "1");
        F3::params.Frobenius_coeffs_c2(1, "5136291436651207728317994048073823738016144056504959939");
        F3::params.Frobenius_coeffs_c2(2, "1073752683758513276629212192812154536507607213288832061");

        F6::params.non_residue("61");
        F6::params.Frobenius_coeffs_c1(0, "1");
        F6::params.Frobenius_coeffs_c1(1, "1073752683758513276629212192812154536507607213288832062");
        F6::params.Frobenius_coeffs_c1(2, "1073752683758513276629212192812154536507607213288832061");
        F6::params.Frobenius_coeffs_c1(3, "6210044120409721004947206240885978274523751269793792000");
        F6::params.Frobenius_coeffs_c1(4, "5136291436651207728317994048073823738016144056504959939");
        F6::params.Frobenius_coeffs_c1(5, "5136291436651207728317994048073823738016144056504959940");
        F2::params.non_residue(F3::params.non_residue()); // F6 cyclotomic_squared()
    }
};

} // namespace snarklib

#endif
