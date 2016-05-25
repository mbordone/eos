/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016 Danny van Dyk
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

#include <eos/form-factors/form-factors.hh>
#include <eos/b-decays/b-to-dstar-l-x-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    template <>
    struct Implementation<BToDstarLeptonInclusiveNeutrinos>
    {
        std::shared_ptr<FormFactors<PToV>> form_factors;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_Dstar;

        UsedParameter m_mu;

        UsedParameter m_tau;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            m_Dstar(p["mass::D^*" + std::string(o.get("q", "d") == "d" ? "+" : "0")], u),
            m_mu(p["mass::mu"], u),
            m_tau(p["mass::tau"], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            if ((o.get("q", "d") != "d") && (o.get("q", "d") != "u")) // q = d is the default
            {
                // only B_{d,u} mesons can decay in this channel
                throw InternalError("BToDstarLeptonInclusiveNeutrinos: q = '" + o["q"] + "' is not a valid option for this decay channel");
            }

            form_factors = FormFactorFactory<PToV>::create("B->D^*@" + o.get("form-factors", "FKN2012"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
        }

        // normalized to N_1 = |V_vb|^2 G_F^2 / (192 pi^3 MB)
        double differential_decay_width_1nu(const double & s) const
        {
            const double mB = m_B(), mB2 = mB * mB;
            const double mDstar = m_Dstar(), mDstar2 = mDstar * mDstar;
            const double lam = lambda(mB2, mDstar2, s), sqrtlam = std::sqrt(lam);

            const double Fperp = std::sqrt(2.0) * sqrtlam / (mB * (mB + mDstar)) * form_factors->v(s);
            const double Fpara = std::sqrt(2.0) * (mB + mDstar) / mB * form_factors->a_1(s);
            const double Flong = 8.0 * mDstar / mB * form_factors->a_12(s);

            return sqrt(lam) * (Flong * Flong * mB2 + (Fpara * Fpara + Fperp * Fperp) * s);
        }

        // normalized to N_1 = |V_vb|^2 G_F^2 / (192 pi^3 MB)
        double normalized_differential_decay_width_1nu(const double & s, const double & z) const
        {
            std::function<double (const double &)> integrand = std::bind(&Implementation<BToDstarLeptonInclusiveNeutrinos>::differential_decay_width_1nu, this, std::placeholders::_1);
            const double s_min = 0.0, s_max = power_of<2>(m_B() - m_Dstar());
            const double Gamma_1 = integrate(integrand, 128, s_min, s_max);

            const double mB = m_B(), mB2 = mB * mB;
            const double mDstar = m_Dstar(), mDstar2 = mDstar * mDstar;
            const double lam = lambda(mB2, mDstar2, s), sqrtlam = std::sqrt(lam);
            const double z2 = z * z;

            const double Fperp = std::sqrt(2.0) * sqrtlam / (mB * (mB + mDstar)) * form_factors->v(s);
            const double Fpara = std::sqrt(2.0) * (mB + mDstar) / mB * form_factors->a_1(s);
            const double Flong = 8.0 * mDstar / mB * form_factors->a_12(s);

            return 3.0 / 8.0 * sqrt(lam) * (
                    2.0 * Flong * Flong * mB2 * (1.0 - z2)
                    + (4.0 * Fperp * Fpara * z + (Fpara * Fpara + Fperp * Fperp) * (1.0 + z2)) * s
                ) / Gamma_1;
        }

        // normalized to N_3 = |V_vb|^2 G_F^2 / (384 pi^3 MB)
        //                   * tau_tau / hbar * G_F^2 m_tau^5 / (192 pi^3)
        double differential_decay_width_3nu(const double & s) const
        {
            const double mB = m_B(), mB2 = mB * mB;
            const double mDstar = m_Dstar(), mDstar2 = mDstar * mDstar;
            const double lam = lambda(mB2, mDstar2, s), sqrtlam = std::sqrt(lam);
            const double mtau = m_tau(), mtau2 = mtau * mtau;
            const double s3 = s * s * s;

            const double Fperp = std::sqrt(2.0) * sqrtlam / (mB * (mB + mDstar)) * form_factors->v(s);
            const double Fpara = std::sqrt(2.0) * (mB + mDstar) / mB * form_factors->a_1(s);
            const double Flong = 8.0 * mDstar / mB * form_factors->a_12(s);
            const double Ftime = sqrtlam / mB2 * form_factors->a_0(s);

            return sqrt(lam) * std::pow(s - mtau2, 2) * (
                    (Fpara * Fpara + Fperp * Fperp) * s * (mtau2 + 2.0 * s)
                    + Flong * Flong * mB2 * (mtau2 + 2.0 * s)
                    + 3.0 * Ftime * Ftime * mB2 * mtau2
                ) / s3;
        }

        // normalized to N_3 = |V_vb|^2 G_F^2 / (384 pi^3 MB)
        //                   * tau_tau / hbar * G_F^2 m_tau^5 / (192 pi^3)
        double normalized_differential_decay_width_3nu(const double & s, const double & snunubar,
                const double & z, const double & phi, const double & zst) const
        {
            std::function<double (const double &)> integrand = std::bind(&Implementation<BToDstarLeptonInclusiveNeutrinos>::differential_decay_width_3nu, this, std::placeholders::_1);
            const double s_min = 3.16, s_max = power_of<2>(m_B() - m_Dstar());
            const double Gamma_3 = integrate(integrand, 128, s_min, s_max);

            const double mB = m_B(), mB2 = mB * mB;
            const double mDstar = m_Dstar(), mDstar2 = mDstar * mDstar;
            const double lam = lambda(mB2, mDstar2, s), sqrtlam = std::sqrt(lam);
            const double z2 = z * z;
            const double mtau = m_tau(), mtau2 = mtau * mtau, mtau4 = mtau2 * mtau2, mtau8 = mtau4 * mtau4;
            const double s2 = s * s, s3 = s2 * s;

            const double Fperp = std::sqrt(2.0) * sqrtlam / (mB * (mB + mDstar)) * form_factors->v(s), Fperp2 = Fperp * Fperp;
            const double Fpara = std::sqrt(2.0) * (mB + mDstar) / mB * form_factors->a_1(s), Fpara2 = Fpara * Fpara;
            const double Flong = 8.0 * mDstar / mB * form_factors->a_12(s), Flong2 = Flong * Flong;
            const double Ftime = sqrtlam / mB2 * form_factors->a_0(s), Ftime2 = Ftime * Ftime;

            // constant in z
            const double a = (mtau2 + 2.0 * snunubar) * ((Fperp2 + Fpara2) * (mtau2 + s) * s + 2.0 * Flong2 * mB2 * s + 2.0 * Ftime2 * mB2 * mtau2)
                           - (mtau2 - 2.0 * snunubar) * ((Fperp2 + Fpara2) * (mtau2 - s) * s - 2.0 * Flong2 * mB2 * s + 2.0 * Ftime2 * mB2 * mtau2) * zst;

            // multiplying z
            const double b = 4.0 * (mtau2 + 2.0 * snunubar) * (Flong * Ftime * mB2 * mtau2 + Fperp * Fpara * s2)
                           - 4.0 * (mtau2 - 2.0 * snunubar) * (Flong * Ftime * mB2 * mtau2 - Fperp * Fpara * s2) * zst;

            // multiplying z^2
            const double c = (mtau2 - s) * (mtau2 + 2.0 * snunubar) * (2.0 * Flong2 * mB2 - (Fperp2 + Fpara2) * s)
                           - (mtau2 + s) * (mtau2 - 2.0 * snunubar) * (2.0 * Flong2 * mB2 - (Fperp2 + Fpara2) * s) * zst;

            // multiplying sqrt(1 - z^2)
            const double d = 4.0 * mtau * sqrt(s) * (mtau2 - 2.0 * snunubar) * (Flong * Ftime * mB2 - Fperp * Fpara * s) * sqrt(1.0 - zst * zst);

            // multiplying z sqrt(1 - z^2)
            const double e = 2.0 * mtau * sqrt(s) * (mtau2 - 2.0 * snunubar) * (2.0 * Flong2 * mB2 - (Fperp2 + Fpara2) * s) * sqrt(1.0 - zst * zst);

            return 3.0 / 8.0 * power_of<2>((mtau2 - s) * (mtau2 - snunubar)) * sqrtlam * (a + b * z + c * z2 + (d + e * z) * sqrt(1.0 - z2) * cos(phi)) / (mtau8 * M_PI * s3) / Gamma_3;
        }
    };

    BToDstarLeptonInclusiveNeutrinos::BToDstarLeptonInclusiveNeutrinos(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDstarLeptonInclusiveNeutrinos>(new Implementation<BToDstarLeptonInclusiveNeutrinos>(parameters, options, *this))
    {
    }

    BToDstarLeptonInclusiveNeutrinos::~BToDstarLeptonInclusiveNeutrinos()
    {
    }

    double
    BToDstarLeptonInclusiveNeutrinos::normalized_differential_decay_width_1nu(const double & s, const double & c_theta_mu) const
    {
        return _imp->normalized_differential_decay_width_1nu(s, c_theta_mu);
    }

    double
    BToDstarLeptonInclusiveNeutrinos::normalized_differential_decay_width_3nu(const double & s, const double & snunubar,
            const double & z, const double & phi, const double & zst) const
    {
        return _imp->normalized_differential_decay_width_3nu(s, snunubar, z, phi, zst);
    }
}
