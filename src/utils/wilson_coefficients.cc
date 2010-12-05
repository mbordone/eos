/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christoph Bobeth
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <src/utils/matrix.hh>
#include <src/utils/model.hh>
#include <src/utils/power_of.hh>
#include <src/utils/qcd.hh>
#include <src/utils/wilson_coefficients.hh>

#include <array>
#include <cmath>
#include <vector>

namespace eos
{
    WilsonCoefficients<BToS> evolve(const std::array<double, 15> & wc_qcd_0,
            const std::array<double, 15> & wc_qcd_1,
            const std::array<double, 15> & wc_qcd_2,
            const double & alpha_s_0, const double & alpha_s, const double & nf, const QCD::BetaFunction & beta)
    {
        using std::array;
        typedef array<array<double, 15>, 15> Matrix;

        // diagonalisation matrix of gamma_qcd_0
        static const Matrix V
        {{
            {{  0, 0, 0, 0, 0.9409878113450298, -0.007502261560628922, 0, 0, 0, 0, -0.00034192561587966793, 0.002764024590670891, -0.8282780998098398, 0, 0 }},
            {{  0, 0, 0, 0, -0.3136626037816766, 0.0025007538535429742, 0, 0, 0, 0, -0.00022795041058644593, 0.0018426830604472593, -0.5521853998732265, 0, 0}},
            {{ -0.00817705933858199, 0, 0, -0.15365781326844577, 0.03485140042018328, 0.09222271788266753,
               -0.8096484594817749, 0.91164065860889, 0, 0,  0.023095095888899648, 0.06294995947461862,
               -0.026294542851110006, -0.725871875481995, 0.0523991277432502}},
            {{ -0.04906235603149206, 0, 0, -0.9416338022084989, -0.1045542012605542, -0.27666815364800085, 0.23516917860724376, -0.3276270369409683, 0, 0, 0.03464264383334918, 0.09442493921192827, -0.03944181427666017, 0.6838800872511366, -0.03929934580743761 }},
            {{  0.0005110662086613764, 0, 0, 0.006938972071662499, -0.008712850105046016, -0.02305567947066682, 0.1355723071850575, -0.12937117286053687, 0, 0, -0.0057737739722248755, -0.015737489868654676, 0.006573635712777292, 0.03686183998099939, -0.0032749454839531324 }},
            {{  0.0030663972519682533, 0, 0, 0.07263235909898472, 0.02613855031513811, 0.06916703841200035, -0.22135905053575103, -0.1496164282505422, 0, 0, -0.00866066095833731, -0.023606234802982035, 0.009860453569165841, -0.049865393718132485, 0.002456209112964845 }},
            {{  0.12265589007872987, 0, 0, 0, 0, -0.27750173826584934, 0, 0, 0, 0, 0.48522496407747906, 0, 0, 0, -0.7859869161487413 }},
            {{  0.7359353404723797, 0, 0, 0, 0, 0.8325052147975449, 0, 0, 0, 0, 0.7278374461162129, 0, 0, 0, 0.5894901871115558 }},
            {{ -0.00766599312992062, 0, 0, 0, 0, 0.06937543456646225, 0, 0, 0, 0, -0.12130624101936928, 0, 0, 0, 0.049124182259296316 }},
            {{ -0.04599595877952373, 0, 0, 0, 0, -0.2081263036993865, 0, 0, 0, 0, -0.18195936152905356, 0, 0, 0, -0.03684313669447223 }},
            {{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.990079849429909, 0, 0, 0 }},
            {{  0, 0, 0, 0, 0, 0, 0, 0, 0.9363291775690445, 1., 0, 0, 0, 0, 0 }},
            {{  0, 0, 0, 0, 0, 0, 0, 0, 0.35112344158839165, 0, 0, 0, 0, 0, 0 }},
            {{  0.6623418064251398, 1., 0, 0.2905020654629326, -0.05702956432394517, -0.3022744864984632, -0.4709407072683506, -0.14983276766672257, 0, 0, 0.43028242387321836, 0.07776353250770997, 0.0816037536758473, 0.03965720499707316, 0.16324343643089229}},
            {{  0, 0, 1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}
        }};

        // inverse of diagonalisation matrix
        static const Matrix V_inverse
        {{
            {{0, 0, 0, 0, 0, 0, 1.2078355801181289, 1.6104474401575037, 4.83134232047253, 6.441789760630009, 0, 0, 0, 0, 0 }},
            {{0.11836177258012, -0.04688112148115692, 0.4728978854883222, 0.140543610478892, 3.57455465561306, -1.5186169423711846, -0.13622538926921474, -1.1495429772208166, 4.200426850366367, -5.351534257227639, -0.10164668578015348, 0, 0, 1, 0 }},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }},
            {{-0.05972704970547934, 0.14183442739855728, -1.1429139456212456, -0.9294558048873448, -5.498186941606234, -0.17448002107637978, -0.18966180496026447, -0.014182080561437292, -2.1820264668086478, 0.75287423482508, 0.06965847045378179, 0, 0, 0, 0 }},
            {{0.7084753475326595, -1.0627130212989893, 0, 0, 0, 0, 0.0031922727161175, -0.0021281818107449983, 0.05107636345787996, -0.03405090897191997, 0, 0, 0, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, 0.40039789229956424, -0.2669319281997093, 6.406366276793023, -4.270910851195349, 0, 0, 0, 0, 0 }},
            {{0.061366169238874266, -0.07345372526830733, 0.2477249556992618, -0.1949894734177815, 5.122876219656581, -2.493259223870642, 0.07527797726126845, -0.0620922336189544, 1.2817327647448542,
              -0.9517042478490128, 0.024794038120004115, 0, 0, 0, 0 }},
            {{-0.014192364080408514, -0.07230371312053796, -0.31773530251832655, -0.21203242136881975, -6.155388277025308, -2.8329981583444037, 0.03666061699520339, -0.0027816034935972023, 0.5151283094744479, -0.007204283660441117, -0.12478967898818555, 0, 0, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.8480012484391772, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1., -2.666666666666667, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, -0.45797771894263023, -0.1526592396475433, -7.32764350308208, -2.4425478343606932, 0, 0, 0, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.010019545974805, 0, 0, 0, 0 }},
            {{-0.402441321833647, -1.2073239655009411, 0, 0, 0, 0, 0.00018906006768086274, 0.00006302002256028741, 0.0030249610828938026, 0.0010083203609645988, 0.0033705090871937156, 0, 0, 0, 0 }},
            {{-0.025035360484996275, -0.04619035883071117, -1.8110811783375629, 0.14795131295326958, -12.280943760100731, -0.7400790922791021, -0.08378645815793977, 0.029891709251547136, -0.22749299097366743, 0.27111400805598945, -0.11165786607759942, 0, 0, 0, 0 }},
            {{0, 0, 0, 0, 0, 0, -1.5078942929387142, 0.25131571548978476, -6.0315771717548765, 1.005262861959134, 0, 0, 0, 0, 0 }}
        }};

        static const Matrix G_qcd_1
        {{
            {{ -257.778, 0., 0., 0., 0., -176.973, 0., 0., 0., 0., 98.4601, 0., 0., 0., -16.0202 }},
            {{ 123.984, -77.3333, 0., 56.1243, 6.05517, 153.89, -20.1579, 105.229, 0., 0., -67.1537, -4.22676, -6.90668, -0.494136, 11.9667 }},
            {{ 0., 0., -77.3333, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }},
            {{ 2.77592, 0., 0., -269.682, -26.5128, -56.7707, 167.181, -186.622, 0., 0., -10.5094, -10.0209, 3.45168, 9.60358, 0.133441 }},
            {{ -0.164451, 0., 0., 0., -8.22222, -1.66533e-16, 0., 0., 0., 0., -0.229923, -0.0195824, 5.86815, 0.,0.060218 }},
            {{ -20.6267, 0., 0., 0., 0., -8.22222, 0., 0., 0., 0., -29.1425, 0., 0., 0., 7.55298 }},
            {{ -4.03975, 0., 0., -4.29724, -4.7255, -10.936, 96.3502, -24.2511, 0., 0., -10.3216, -3.3874, 2.04999, 14.889, 1.58061 }},
            {{ 1.41375, 0., 0., 18.2442, 22.9293, 64.6564, -215.241, 74.8239, 0., 0., -2.18519, -10.4539, -0.462874, -43.6055, 4.4996 }},
            {{ 1.25757, 0., 0., 19.9161, -10.3227, -54.7136, 170.721, -120.194, 73.1481, 0., 6.31785, 17.183, 1.19819, 10.437, -4.36506 }},
            {{ -4.80811, 0., 0., -21.958, 6.55206, 34.728, -161.111, 130.879, 12.1723, 96.2963, -23.5699, -23.6213, -4.47007, -11.5095, 3.22255 }},
            {{ -17.8632, 0., 0., 0., 0., -156.32, 0., 0., 0., 0., -40.5556, 0., 0., 0., -28.0772 }},
            {{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 36.1111, 0., 0., 0. }},
            {{ 0.00737421, 0., 0., 0., 31.0528, -0.183045, 0., 0., 0., 0., -1.38778e-17, 0.255842, -40.5556, 0., 0.0115907 }},
            {{ 2.2264, 0., 0., 35.1232, -4.69861, 0.628604, 190.083, -289.392, 0., 0., -16.0293, -27.9011, 8.65342, 77.432, 1.76299 }},
            {{ -3.12107, 0., 0., 0., 0., 190.81, 0., 0., 0., 0., -175.6, 0., 0., 0., 108.722 }}
        }};

        static const double zeta_3 = 1.2020569031595943;
        double u11 = -1927.0 / 2 + 257.0 / 9 * nf + 40.0 / 9 * nf * nf + (224 + 160.0 / 3 * nf) * zeta_3;
        double u12 = 475.0 / 9 + 362.0 / 27 * nf - 40.0 / 27 * nf * nf - (896.0 / 3 + 320.0 / 9 * nf) * zeta_3;
        double u21 = 307.0 / 2 + 361.0 / 3 * nf - 20.0 / 3 * nf * nf - (1344 + 160* nf) * zeta_3;
        double u22 = 1298.0 / 3 - 76.0 / 3 * nf - 224* zeta_3;
        double u13 = 269107.0 / 13122 - 2288.0 / 729 * nf - 1360.0 / 81* zeta_3;
        double u14 = -2425817.0 / 13122 + 30815.0 / 4374 * nf - 776.0 / 81* zeta_3;
        double u23 = 69797.0 / 2187 + 904.0 / 243 * nf + 2720.0 / 27* zeta_3;
        double u24 = 1457549.0 / 8748 - 22067.0 / 729 * nf - 2768.0 / 27* zeta_3;
        double u33 = -4203068.0 / 2187 + 14012.0 / 243 * nf - 608.0 / 27* zeta_3;
        double u34 = -18422762.0 / 2187 + 888605.0 / 2916 * nf + 272.0 / 27 * nf * nf
                    + (39824.0 / 27 + 160. * nf) * zeta_3;
        double u43 = -5875184.0 / 6561 + 217892.0 / 2187 * nf + 472.0 / 81 * nf * nf
                    + (27520.0 / 81 + 1360.0 / 9 * nf) * zeta_3;
        double u44 = -70274587.0 / 13122 + 8860733.0 / 17496 * nf - 4010.0 / 729 * nf * nf
                    + (16592.0 / 81 + 2512.0 / 27 * nf) * zeta_3;
        double u53 = -194951552.0 / 2187 + 358672.0 / 81 * nf - 2144.0 / 81 * nf * nf + 87040.0 / 27* zeta_3;
        double u54 = -130500332.0 / 2187 - 2949616.0 / 729 * nf + 3088.0 / 27 * nf * nf
                    + (238016.0 / 27 + 640. * nf) * zeta_3;
        double u63 = 162733912.0 / 6561 - 2535466.0 / 2187 * nf + 17920.0 / 243 * nf * nf
                    + (174208.0 / 81 + 12160.0 / 9 * nf) * zeta_3;
        double u64 = 13286236.0 / 6561 - 1826023.0 / 4374 * nf - 159548.0 / 729 * nf * nf
                    - (24832.0 / 81 + 9440.0 / 27 * nf) * zeta_3;
        double u15 = -343783.0 / 52488 + 392.0 / 729 * nf + 124.0 / 81* zeta_3;
        double u16 = -37573.0 / 69984 + 35.0 / 972 * nf + 100.0 / 27* zeta_3;
        double u25 = -37889.0 / 8748 - 28.0 / 243 * nf - 248.0 / 27* zeta_3;
        double u26 = 366919.0 / 11664 - 35.0 / 162 * nf - 110.0 / 9* zeta_3;
        double u35 = 674281.0 / 4374 - 1352.0 / 243 * nf - 496.0 / 27* zeta_3;
        double u36 = 9284531.0 / 11664 - 2798.0 / 81 * nf - 26.0 / 27* nf * nf
                    - (1921.0 / 9 + 20* nf) * zeta_3;
        double u45 = 2951809.0 / 52488 - 31175.0 / 8748 * nf - 52.0 / 81* nf * nf
                    - (3154.0 / 81 + 136.0 / 9* nf) * zeta_3;
        double u46 = 3227801.0 / 8748 - 105293.0 / 11664 * nf - 65.0 / 54* nf * nf
                    + (200.0 / 27 - 220.0 / 9* nf) * zeta_3;
        double u55 = 14732222.0 / 2187 - 27428.0 / 81 * nf + 272.0 / 81* nf * nf
                    - 13984.0 / 27* zeta_3;
        double u56 = 16521659.0 / 2916 + 8081.0 / 54 * nf - 316.0 / 27* nf * nf
                    - (22420.0 / 9 + 200* nf) * zeta_3;
        double u65 = -22191107.0 / 13122 + 395783.0 / 4374 * nf - 1720.0 / 243* nf * nf
                    - (33832.0 / 81 + 1360.0 / 9 * nf) * zeta_3;
        double u66 = -32043361.0 / 8748 + 3353393.0 / 5832 * nf - 533.0 / 81* nf * nf
                    + (9248.0 / 27 - 1120.0 / 9* nf) * zeta_3;
        static const double u17 = -13234.0 / 2187;
        static const double u18 = 13957.0 / 2916;
        static const double u19 = -1359190.0 / 19683 + 6976.0 / 243 * zeta_3;
        static const double u27 = 20204.0 / 729;
        static const double u28 = 14881.0 / 972;
        static const double u29 = -229696.0 / 6561 - 3584.0 / 81 * zeta_3;
        static const double u37 = 92224.0 / 729;
        static const double u38 = 66068.0 / 243;
        static const double u39 = -1290092.0 / 6561 + 3200.0 / 81 * zeta_3;
        static const double u47 = -184190.0 / 2187;
        static const double u48 = -1417901.0 / 5832;
        static const double u49 = -819971.0 / 19683 - 19936.0 / 243 * zeta_3;
        static const double u57 = 1571264.0 / 729;
        static const double u58 = 3076372.0 / 243;
        static const double u59 = -16821944.0 / 6561 + 30464.0 / 81 * zeta_3;
        static const double u67 = -1792768.0 / 2187;
        static const double u68 = -3029846.0 / 729;
        static const double u69 = -17787368.0 / 19683 - 286720.0 / 243 * zeta_3;
        static const double u99 = -9769.0 / 27;
        Matrix gamma_qcd_2_transposed
        {{
            {{ u11, u21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u12, u22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u13, u23, u33, u43, u53, u63, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u14, u24, u34, u44, u54, u64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u15, u25, u35, u45, u55, u65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u16, u26, u36, u46, u56, u66, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u17, u27, u37, u47, u57, u67, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u18, u28, u38, u48, u58, u68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }},
            {{ u19, u29, u39, u49, u59, u69, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, u99, 0.0 }},
            {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, u99 }},
        }};
        Matrix G_qcd_2 = V_inverse * gamma_qcd_2_transposed * V;

        static const array<double, 15> gamma_qcd_0_eigenvalues
        {{
            -16.000000000000000, -15.3333333333333334, -15.333333333333334, -13.790720905057988, -8.000000000000005,
            - 7.999999999999997, - 6.4858257980688645, + 6.265491004149404, - 6.000000000000000, -4.666666666666667,
            + 4.000000000000004, + 4.0000000000000000, + 3.999999999999999, + 2.233277921199660, +2.000000000000003
        }};

        WilsonCoefficients<BToS> result;
        result._alpha_s = alpha_s;
        double eta = alpha_s_0 / alpha_s;

        array<double, 15> a;
        Matrix H_qcd_0, H_qcd_1, H_qcd_2;
        for (unsigned i(0) ; i < a.size() ; ++i)
        {
            a[i] = gamma_qcd_0_eigenvalues[i] / 2.0 / beta[0];
        }

        for (unsigned i(0) ; i < a.size() ; ++i)
        {
            for (unsigned j(0) ; j < a.size() ; ++j)
            {
                H_qcd_0[i][j] = 0.0;
                H_qcd_1[i][j] = -G_qcd_1[i][j] / (2.0 * beta[0]) / (1.0 + a[i] - a[j]);
            }
            H_qcd_0[i][i] = std::pow(eta, a[i]);
            H_qcd_1[i][i] += beta[1] / beta[0] * a[i];
        }

        // Need complete H_qcd_1 to compute H_qcd_2!
        for (unsigned i(0) ; i < a.size() ; ++i)
        {
            for (unsigned j(0) ; j < a.size() ; ++j)
            {
                H_qcd_2[i][j] = -G_qcd_2[i][j] / (2.0 * beta[0]) / (2.0 + a[i] - a[j]);
                H_qcd_2[i][j] += -beta[1] / beta[0] * (1.0 + a[i] - a[j]) / (2.0 + a[i] - a[j]) * H_qcd_1[i][j];

                for (unsigned k(0) ; k < a.size() ; ++k)
                {
                    H_qcd_2[i][j] += (1.0 + a[i] - a[k]) / (2.0 + a[i] - a[j]) * H_qcd_1[i][k] * H_qcd_1[k][j];
                }
            }

            H_qcd_2[i][i] += beta[2] / 2.0 / beta[0] * a[i];
        }

        Matrix U_qcd_0 = V * H_qcd_0 * V_inverse;
        Matrix U_qcd_1 = V * ((H_qcd_1 * H_qcd_0) + (-eta) * (H_qcd_0 * H_qcd_1)) * V_inverse;
        Matrix U_qcd_2 = V * (H_qcd_2 * H_qcd_0 + (-eta) * H_qcd_1 * H_qcd_0 * H_qcd_1 + (-eta*eta) * H_qcd_0 * (H_qcd_2 - H_qcd_1 * H_qcd_1)) * V_inverse;
        double a_s = alpha_s / (4 * M_PI);
        array<double, 15> result_qcd_0 = U_qcd_0 * wc_qcd_0;
        array<double, 15> result_qcd_1 = U_qcd_1 * wc_qcd_0 + eta * (U_qcd_0 * wc_qcd_1);
        array<double, 15> result_qcd_2 = U_qcd_2 * wc_qcd_0 + eta * (U_qcd_1 * wc_qcd_1) + power_of<2>(eta) * (U_qcd_0 * wc_qcd_2);

        result._coefficients = result_qcd_0 + a_s * result_qcd_1 + power_of<2>(a_s) * result_qcd_2;

        return result;
    }
}
