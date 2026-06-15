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
    void define_fracture_network_description( pybind11::module& module )
    {
        pybind11::class_< FractureSamplerConfig >( module,
            "FractureSamplerConfig",
            "Configuration of the stochastic fracture sampler." )
            .def( pybind11::init<>() )
            .def_readwrite( "move_ratio", &FractureSamplerConfig::move_ratio,
                "Relative probability of proposing a change move." )
            .def_readwrite( "length", &FractureSamplerConfig::length,
                "Distribution used to sample fracture lengths." )
            .def_readwrite( "azimuth", &FractureSamplerConfig::azimuth,
                "Distribution used to sample fracture orientations." );

        pybind11::class_< FractureSetDescription >( module,
            "FractureSetDescription", "Description of a fracture family." )
            .def( pybind11::init<>() )

            .def_readwrite( "name", &FractureSetDescription::fset_name,
                "Name of the fracture set." )

            .def_readwrite( "sampler", &FractureSetDescription::sampler,
                "Sampling configuration used to generate fractures." )

            .def_readwrite( "birth_ratio", &FractureSetDescription::birth_ratio,
                "Relative probability of birth moves." )

            .def_readwrite( "death_ratio", &FractureSetDescription::death_ratio,
                "Relative probability of death moves." )

            .def_readwrite( "change_ratio",
                &FractureSetDescription::change_ratio,
                "Relative probability of change moves." )

            .def_readwrite( "p20", &FractureSetDescription::p20,
                "Target fracture density (number of fractures per unit area)." )

            .def_readwrite( "expected_number",
                &FractureSetDescription::expected_number,
                "Expected number of fractures, if monitored." )

            .def_readwrite( "p21", &FractureSetDescription::p21,
                "Target fracture intensity (total fracture length per unit "
                "area)." )

            .def_readwrite( "expected_total_length",
                &FractureSetDescription::expected_total_length,
                "Expected cumulative fracture length, if monitored." )

            .def(
                "add_observed_fracture",
                []( FractureSetDescription& self, const geode::Point2D& start,
                    const geode::Point2D& end ) {
                    self.observed_fractures.push_back( { start, end } );
                },
                "Add a fixed observed fracture defined by two endpoints." )

            .def_readwrite( "minimal_spacing",
                &FractureSetDescription::minimal_spacing,
                "Minimum allowed spacing between fractures." )

            .def( "string", &FractureSetDescription::string,
                "Return a detailed textual description of the fracture set." )

            .def( "__repr__", []( const FractureSetDescription& self ) {
                return "<FractureSetDescription name='" + self.fset_name + "'>";
            } );

        pybind11::class_< FractureNetworkDescription >( module,
            "FractureNetworkDescription",
            "Description of a complete fracture network simulation." )

            .def( pybind11::init<>() )

            .def_readwrite( "name", &FractureNetworkDescription::fnet_name,
                "Name of the fracture network." )

            .def_readwrite( "domain", &FractureNetworkDescription::domain,
                "Spatial simulation domain." )

            .def_readwrite( "fracture_sets",
                &FractureNetworkDescription::fracture_sets,
                "List of fracture sets composing the network." )

            .def( "add_fracture_set",
                &FractureNetworkDescription::add_fracture_set,
                pybind11::return_value_policy::reference_internal,
                "Create and return a new fracture set." )

            .def( "add_x_node_monitoring",
                &FractureNetworkDescription::add_x_node_monitoring,
                pybind11::arg( "beta" ),
                "Enable monitoring of X-node interactions with the given "
                "weight." )

            .def_readwrite( "expected_x_node",
                &FractureNetworkDescription::expected_x_node,
                "Expected number of X-nodes, if monitored." )

            .def( "string", &FractureNetworkDescription::string,
                "Return a detailed textual description of the network." )

            .def( "__repr__", []( const FractureNetworkDescription& self ) {
                return "<FractureNetworkDescription name='" + self.fnet_name
                       + "'>";
            } );

        //        module.def( "build_fractures_simulation_context",
        //            &build_fractures_simulation_context, pybind11::arg(
        //            "description" ), "Build a simulation context from a
        //            fracture network description." );
        //
        //        module.def( "build_fractures_targeted_stat",
        //            &build_fractures_targeted_stat, pybind11::arg(
        //            "description" ), "Create the target statistics associated
        //            with a fracture network " "description." );

        module.def( "build_fractures_simulation_runner",
            &build_fractures_simulation_runner, pybind11::arg( "description" ),
            "Create a ready-to-use simulation runner from a fracture network "
            "description." );
    }

} // namespace geode