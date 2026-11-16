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
            "Definition of an object set of a simulation: its name, the "
            "sampling of its objects, the relative probabilities of its "
            "moves and its fixed objects. Returned by the builder of the "
            "simulation, which it keeps alive." );
        pyclass.def_readonly( "name", &Definition::name, "Object set name." );
        def_readwrite_named( pyclass, "sampler", &Definition::sampler,
            "Sampling configuration used to generate objects." );
        def_readwrite_named( pyclass, "dynamics", &Definition::dynamics,
            "Relative probabilities of the birth / death / change moves." );
        pyclass
            .def(
                "nb_fixed_objects",
                []( const Definition& self ) {
                    return self.fixed_objects.size();
                },
                "Return the number of fixed (observed) objects." )
            .def( "__repr__", [pyclass_name]( const Definition& self ) {
                return python_call_repr( pyclass_name,
                    { "name=" + python_repr( self.name ),
                        "nb_fixed_objects="
                            + std::to_string( self.fixed_objects.size() ) } );
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
                "Add a fixed (observed) segment from start to end: the "
                "simulation never moves nor removes it, and it comes first "
                "in the segments of the set." );
    }
} // namespace geode
