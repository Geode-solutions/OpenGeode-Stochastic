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
            "Builder of a fracture network simulation (density, intensity, "
            "spacing and X-node terms)." )
            .def( pybind11::init<>() )
            .def( "set_domain", &FractureProcessBuilder::set_domain,
                pybind11::arg( "domain" ),
                "Set the spatial simulation domain." )
            .def( "add_fracture_set", &FractureProcessBuilder::add_fracture_set,
                pybind11::arg( "name" ), pybind11::arg( "p20" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                pybind11::return_value_policy::reference_internal,
                "Add a fracture set with a density term (number of fractures "
                "per unit area) and return its definition." )
            .def( "add_intensity", &FractureProcessBuilder::add_intensity,
                pybind11::arg( "set_name" ), pybind11::arg( "p21" ),
                pybind11::arg( "expected_total_length" ) = std::nullopt,
                "Add an intensity term (fracture length per unit area) on a "
                "fracture set." )
            .def( "add_minimal_spacing",
                &FractureProcessBuilder::add_minimal_spacing,
                pybind11::arg( "set_name" ), pybind11::arg( "minimal_spacing" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                "Forbid fractures of a set closer than the given distance (0 "
                "forbids intersections)." )
            .def( "add_x_node_interaction",
                &FractureProcessBuilder::add_x_node_interaction,
                pybind11::arg( "set_names" ), pybind11::arg( "beta" ),
                pybind11::arg( "expected_count" ) = std::nullopt,
                "Penalize intersections between fractures of different sets "
                "(beta in [0, 1])." )
            .def( "expected_statistics",
                &FractureProcessBuilder::expected_statistics,
                pybind11::return_value_policy::copy,
                "Expected statistics registered with the model terms." )
            .def(
                "build_simulation_runner",
                []( const FractureProcessBuilder& self ) {
                    return FractureSimulationRunner{
                        self.build_simulation_context()
                    };
                },
                "Create a ready-to-use simulation runner." );
    }
} // namespace geode
