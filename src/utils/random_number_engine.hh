/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH
#define WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH 1

#include <src/utils/private_implementation_pattern.hh>

#include <cstdint>

namespace eos
{
    /*!
     * Generate pseudo-random floating point numbers
     * in the range 0 ... 1.
     *
     * RandomNumberEngine keeps its state/seed on copying.
     */
    class RandomNumberEngine :
        public PrivateImplementationPattern<RandomNumberEngine>
    {
        public:
            ///@{
            /// Basic functions.
            RandomNumberEngine();
            ~RandomNumberEngine();
            ///@}

            /// Obtain a pseudo-random number in the range 0.0 .. 1.0.
            double operator() ();

            /// Return the maximal value that can be drawn.
            double max() const { return 1.0; }

            /// Return the minimal value that can be drawn.
            double min() const { return 0.0; }
    };
}

#endif
