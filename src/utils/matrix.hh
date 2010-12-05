/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#ifndef EOS_GUARD_SRC_UTILS_MATRIX_HH
#define EOS_GUARD_SRC_UTILS_MATRIX_HH 1

#include <array>

namespace eos
{
    /* Addition */

    /* matrix plus matrix */
    template <std::size_t m_, std::size_t n_>
    std::array<std::array<double, n_>, m_> operator+ (const std::array<std::array<double, n_>, m_> & x,
            const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] + y[i][j];
            }
        }

        return result;
    }

    /* matrix minus matrix */
    template <std::size_t m_, std::size_t n_>
    std::array<std::array<double, n_>, m_> operator- (const std::array<std::array<double, n_>, m_> & x,
            const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] - y[i][j];
            }
        }

        return result;
    }

    /* vector plus vector */
    template <std::size_t m_>
    std::array<double, m_> operator+ (const std::array<double, m_> & x,
            const std::array<double, m_> & y)
    {
        std::array<double, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            result[i] = x[i] + y[i];
        }

        return result;
    }

    /* Multiplication */

    /* matrix times matrix */
    template <std::size_t m_, std::size_t n_, std::size_t o_>
    std::array<std::array<double, n_>, m_> operator* (const std::array<std::array<double, o_>, m_> & x,
            const std::array<std::array<double, n_>, o_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = 0.0;

                for (std::size_t k(0) ; k < o_ ; ++k)
                {
                    result[i][j] += x[i][k] * y[k][j];
                }
            }
        }

        return result;
    }

    /* matrix times vector */
    template <std::size_t m_, std::size_t n_>
    std::array<double, m_> operator* (const std::array<std::array<double, n_>, m_> & x, const std::array<double, n_> & y)
    {
        std::array<double, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            result[i] = 0.0;

            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i] += x[i][j] * y[j];
            }
        }

        return result;
    }

    /* scalar times matrix */
    template <std::size_t m_, std::size_t n_>
    std::array<std::array<double, n_>, m_> operator* (const double & x, const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result = y;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] *= x;
            }
        }

        return result;
    }

    /* scalar times vector */
    template <std::size_t n_>
    std::array<double, n_> operator* (const double & x, const std::array<double, n_> & y)
    {
        std::array<double, n_> result = y;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result[i] *= x;
        }

        return result;
    }
}

#endif
