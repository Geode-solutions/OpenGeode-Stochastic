/*
 * Copyright (c) 2019 - 2026 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#pragma once

#include "../../../../common.hpp"

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>

namespace geode
{
    template < typename ObjectType >
    pybind11::class_< ObjectSetDefinition< ObjectType > >
        define_object_set_definition(
            pybind11::module& module, const std::string& typestr )
    {
        using Definition = ObjectSetDefinition< ObjectType >;
        const auto pyclass_name = absl::StrCat( typestr, "SetDefinition" );

        pybind11::class_< Definition > pyclass( module, pyclass_name.c_str(),
            "Sampling definition of an object set." );
        pyclass.def_readonly( "name", &Definition::name, "Object set name." )
            .def_readwrite( "sampler", &Definition::sampler,
                "Sampling configuration used to generate objects." )
            .def_readwrite( "dynamics", &Definition::dynamics,
                "Relative probabilities of the birth / death / change moves." )
            .def(
                "nb_fixed_objects",
                []( const Definition& self ) {
                    return self.fixed_objects.size();
                },
                "Number of fixed (observed) objects." )
            .def( "__repr__", [pyclass_name]( const Definition& self ) {
                return absl::StrCat(
                    "<", pyclass_name, " name='", self.name, "'>" );
            } );
        return pyclass;
    }

    void define_simulation_context( pybind11::module& module )
    {
        define_object_set_definition< OwnerSegment2D >( module, "Segment2D" )
            .def(
                "add_fixed_segment",
                []( ObjectSetDefinition< OwnerSegment2D >& self,
                    const Point2D& start, const Point2D& end ) {
                    self.fixed_objects.emplace_back( start, end );
                },
                pybind11::arg( "start" ), pybind11::arg( "end" ),
                "Add a fixed (observed) segment defined by two endpoints." );
    }
} // namespace geode
