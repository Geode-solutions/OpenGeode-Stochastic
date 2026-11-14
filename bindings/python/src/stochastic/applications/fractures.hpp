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

#include "../../common.hpp"
#include <geode/stochastic/applications/fractures.hpp>

namespace geode
{
    void define_fracture_process_builder( pybind11::module& module )
    {
        pybind11::class_< FractureProcessBuilder >( module,
            "FractureProcessBuilder",
            "Builder of a fracture network simulation: a 2D domain, fracture "
            "sets (segments) and terms of the model (density, intensity, "
            "minimal spacing and X-node interactions). Each term has a "
            "statistic, whose expected mean can be given to validate the "
            "simulation (see expected_statistics)." )
            .def( pybind11::init<>(), "Create a builder with no domain and no "
                                      "fracture set." )
            .def( "__repr__",
                []( const FractureProcessBuilder& ) {
                    return "FractureProcessBuilder()";
                } )
            .def( "set_domain", &FractureProcessBuilder::set_domain,
                pybind11::arg( "domain" ),
                "Set the spatial domain (a SpatialDomainConfig2D, copied); it "
                "is checked when the simulation is built." )
            .def( "add_fracture_set", &FractureProcessBuilder::add_fracture_set,
                pybind11::arg( "name" ), pybind11::arg( "p20" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                pybind11::return_value_policy::reference_internal,
                "Add a fracture set with a density term named name + "
                "\"_p20\" (p20: number of fractures per unit area) and "
                "return its Segment2DSetDefinition, in which the sampling "
                "and the fixed fractures are set. expected_count, if given, "
                "is the expected mean of the statistic of the term (the "
                "number of fractures in the domain)." )
            .def( "add_intensity", &FractureProcessBuilder::add_intensity,
                pybind11::arg( "set_name" ), pybind11::arg( "p21" ),
                pybind11::arg( "expected_total_length" ) = std::nullopt,
                "Add an intensity term named set_name + \"_p21\" (p21: "
                "fracture length per unit area) on a fracture set. "
                "expected_total_length, if given, is the expected mean of "
                "the statistic of the term (the total length of the "
                "fractures inside the domain). The set is looked up when the "
                "simulation is built." )
            .def( "add_minimal_spacing",
                &FractureProcessBuilder::add_minimal_spacing,
                pybind11::arg( "set_name" ), pybind11::arg( "minimal_spacing" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                "Add a hard-core term named set_name + \"_spacing\": "
                "fractures of the set cannot be closer than minimal_spacing "
                "(0 forbids intersections). expected_count, if given, is the "
                "expected mean of the statistic of the term (the number of "
                "interacting pairs in the domain). The set is looked up when "
                "the simulation is built." )
            .def( "add_x_node_interaction",
                &FractureProcessBuilder::add_x_node_interaction,
                pybind11::arg( "set_names" ), pybind11::arg( "beta" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                "Add a term penalizing the intersections (X-nodes) between "
                "fractures of different sets among set_names: beta in [0, "
                "1], 0 forbids X-nodes and 1 has no effect. Raise an "
                "OpenGeodeDataError if beta is out of [0, 1] or fewer than "
                "two sets are given. expected_count, if given, is the "
                "expected mean of the statistic of the term (the number of "
                "intersections in the domain)." )
            .def( "expected_statistics",
                &FractureProcessBuilder::expected_statistics,
                pybind11::return_value_policy::copy,
                "Return the expected statistics given with the terms, as a "
                "list of TargetStatisticConfig, to give to "
                "validate_statistics." )
            .def(
                "build_simulation_runner",
                []( const FractureProcessBuilder& self ) {
                    return FractureSimulationRunner{
                        self.build_simulation_context()
                    };
                },
                pybind11::call_guard< pybind11::gil_scoped_release >(),
                "Create a simulation runner from the current definition (the "
                "builder can still be modified, which does not change the "
                "runner). Raise an OpenGeodeDataError if the domain, a set or "
                "a term is invalid." );
    }
} // namespace geode
