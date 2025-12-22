/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/stochastic/sampling/mcmc/helpers/fracture_simulation_runner.hpp>

namespace geode
{
    void define_fracture_simulation( pybind11::module& module )
    {
        using namespace geode;

        pybind11::class_< FractureSetDescription >(
            module, "FractureSetDescription" )
            .def( pybind11::init<>() )
            .def_readwrite( "name", &FractureSetDescription::name )
            .def_readwrite( "length", &FractureSetDescription::length )
            .def_readwrite( "azimuth", &FractureSetDescription::azimuth )
            .def_readwrite( "p20", &FractureSetDescription::p20 )
            .def_readwrite( "p20", &FractureSetDescription::p21 )
            .def( "add_observed_fracture",
                []( FractureSetDescription& self, const geode::Point2D& a,
                    const geode::Point2D& b ) {
                    self.observed_fractures.push_back( { a, b } );
                } )
            .def_readwrite(
                "minimal_spacing", &FractureSetDescription::minimal_spacing )
            .def_readwrite(
                "birth_ratio", &FractureSetDescription::birth_ratio )
            .def_readwrite(
                "death_ratio", &FractureSetDescription::death_ratio )
            .def_readwrite(
                "change_ratio", &FractureSetDescription::change_ratio )
            .def( "string", &FractureSetDescription::string,
                "Return a detailed description of the fracture set" )
            .def( "__repr__", []( const FractureSetDescription& self ) {
                return "<FractureSetDescription name='" + self.name + "'>";
            } );

        pybind11::class_< FractureSimulationRunner,
            std::shared_ptr< FractureSimulationRunner > >(
            module, "FractureSimulationRunner" )
            .def( pybind11::init< const SpatialDomain< 2 >& >(),
                pybind11::arg( "box" ) )
            .def( "add_x_node_monitoring",
                &FractureSimulationRunner::add_x_node_monitoring,
                pybind11::arg( "double" ),
                "Add a monitoring value for x node, value should be "
                "in[0.,1.]." )
            .def( "add_fracture_set_descriptor",
                &FractureSimulationRunner::add_fracture_set_descriptor,
                pybind11::arg( "descriptor" ),
                "Add a fracture set configuration to the simulation." )
            .def( "initialize", &FractureSimulationRunner::initialize,
                "Initialize internal samplers, energy terms, and proposal "
                "kernels." )
            .def( "check_statistics",
                &FractureSimulationRunner::check_statistics,
                pybind11::arg( "statistic_monitoring" ),
                "Check computed statistics after simulation." )
            .def( "string", &FractureSimulationRunner::string,
                "Return a detailed description of the simulation "
                "configurator." )
            // Explicit overload bindings
            //            .def( "run",
            //                static_cast< const ObjectSets< OwnerSegment2D >& (
            //                    FractureSimulationRunner::*) ( RandomEngine&,
            //                    index_t ) >( &FractureSimulationRunner::run ),
            //                pybind11::arg( "engine" ), pybind11::arg( "steps"
            //                ), "Run simulation for a fixed number of steps." )
            .def( "run",
                static_cast< StatisticsMonitor ( FractureSimulationRunner::* )(
                    RandomEngine&, const SimulationConfigurator& ) >(
                    &FractureSimulationRunner::run ),
                pybind11::arg( "engine" ), pybind11::arg( "config" ),
                "Run the simulation and return statistics monitoring "
                "results." );
    }
} // namespace geode