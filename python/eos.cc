/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#include "eos/observable.hh"
#include "eos/utils/kinematic.hh"
#include "eos/utils/parameters.hh"
#include "eos/utils/options.hh"

#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

using namespace boost::python;
using namespace eos;

namespace impl
{
    // raw constructor for class Kinematics
    object
    Kinematics_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        list items = kwargs.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
            object value = items[i][1];
            self.attr("declare")(name, value);
        }

        return object();
    }

    // raw constructor for class Options
    object
    Options_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        list items = kwargs.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
            object value = items[i][1];
            self.attr("set")(name, value);
        }

        return object();
    }
}

BOOST_PYTHON_MODULE(eos)
{
    using namespace boost::python;
    using namespace eos;

    // Parameters
    class_<Parameters>("Parameters", no_init)
        .def("Defaults", &Parameters::Defaults)
        .staticmethod("Defaults")
        .def("__getitem__", (Parameter (Parameters::*)(const std::string &) const) &Parameters::operator[])
        ;

    // Parameter
    class_<Parameter>("Parameter", no_init)
        .def(float_(self))
        .def("name", &Parameter::name, return_value_policy<copy_const_reference>())
        .def("set", &Parameter::set)
        ;

    // Kinematics
    class_<Kinematics>("Kinematics", no_init)
        .def("__init__", raw_function(&impl::Kinematics_ctor))
        .def(init<>())
        .def("__getitem__", (KinematicVariable (Kinematics::*)(const std::string &) const) &Kinematics::operator[])
        .def("declare", (KinematicVariable (Kinematics::*)(const std::string &, const double &)) &Kinematics::declare, return_value_policy<return_by_value>())
        .def("as_string", &Kinematics::as_string)
        ;

    // KinematicVariable
    class_<KinematicVariable>("KinematicVariable", no_init)
        .def(float_(self))
        .def("name", &KinematicVariable::name, return_value_policy<copy_const_reference>())
        .def("set", &KinematicVariable::set)
        ;

    // Options
    class_<Options>("Options", no_init)
        .def("__init__", raw_function(&impl::Options_ctor))
        .def(init<>())
        .def("set", &Options::set)
        .def("as_string", &Options::as_string)
        ;

    register_ptr_to_python<std::shared_ptr<Observable>>();

    // Observable
    class_<Observable, boost::noncopyable>("Observable", no_init)
        .def("make", &Observable::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        .def("evaluate", &Observable::evaluate)
        .def("name", &Observable::name, return_value_policy<copy_const_reference>())
        ;
}
