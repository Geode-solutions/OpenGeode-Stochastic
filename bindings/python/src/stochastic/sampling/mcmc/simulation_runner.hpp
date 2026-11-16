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

#include "../../../common.hpp"

#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>

#include <geode/geometry/basic_objects/segment.hpp>

namespace geode
{
    void define_simulation_configurator( pybind11::module& module )
    {
        using SimulationConfigurator = geode::SimulationConfigurator;
        pybind11::class_< SimulationConfigurator > config( module,
            "SimulationConfigurator",
            "Parameters of a simulation: after the burn-in steps, each of "
            "the realizations is the state reached after the given number of "
            "Metropolis-Hastings steps from the previous one." );
        config.def(
            pybind11::init(
                []( std::optional< index_t > realizations,
                    std::optional< index_t > metropolis_hasting_steps,
                    std::optional< index_t > burn_in_steps,
                    std::optional< SimulationPrinterConfigurator > printer ) {
                    SimulationConfigurator config;
                    detail::set_if_given( config.realizations, realizations );
                    detail::set_if_given( config.metropolis_hasting_steps,
                        metropolis_hasting_steps );
                    detail::set_if_given( config.burn_in_steps, burn_in_steps );
                    config.printer = std::move( printer );
                    return config;
                } ),
            pybind11::kw_only(), pybind11::arg( "realizations" ) = std::nullopt,
            pybind11::arg( "metropolis_hasting_steps" ) = std::nullopt,
            pybind11::arg( "burn_in_steps" ) = std::nullopt,
            pybind11::arg( "printer" ) = std::nullopt,
            "Create a configuration (defaults: 1000 realizations, 1000 "
            "steps between realizations, 1000 burn-in steps, no "
            "printer)." );
        def_readwrite_named( config, "realizations",
            &SimulationConfigurator::realizations,
            "Number of realizations to generate." );
        def_readwrite_named( config, "metropolis_hasting_steps",
            &SimulationConfigurator::metropolis_hasting_steps,
            "Number of Metropolis-Hastings steps before each realization." );
        def_readwrite_named( config, "burn_in_steps",
            &SimulationConfigurator::burn_in_steps,
            "Number of burn-in steps before the first realization." );
        def_readwrite_named( config, "printer",
            &SimulationConfigurator::printer,
            "SimulationPrinterConfigurator of the output files, None (the "
            "default) to write nothing. To stop printing, create a "
            "SimulationConfigurator without printer." );
        config
            .def( "string", &SimulationConfigurator::string,
                "Return a readable description of the configuration." )
            .def( "__repr__", []( const SimulationConfigurator& self ) {
                return python_call_repr( "SimulationConfigurator",
                    { "realizations=" + python_repr( self.realizations ),
                        "metropolis_hasting_steps="
                            + python_repr( self.metropolis_hasting_steps ),
                        "burn_in_steps=" + python_repr( self.burn_in_steps ),
                        "printer=" + python_repr( self.printer ) } );
            } );
    }

    template < typename object_type >
    void define_simulation_runner(
        pybind11::module_& module, const std::string& typestr )
    {
        using Runner = geode::SimulationRunner< object_type >;
        const auto pyclass_name = absl::StrCat( typestr, "SimulationRunner" );

        pybind11::class_< Runner >( module, pyclass_name.c_str(),
            "Metropolis-Hastings simulation of object sets, created by a "
            "builder (see FractureProcessBuilder.build_simulation_runner). "
            "The calls running the simulation release the GIL: do not use "
            "the runner and its results, nor the RandomEngine given to the "
            "call, from another thread meanwhile." )
            .def( "run",
                pybind11::overload_cast< geode::RandomEngine&, geode::index_t >(
                    &Runner::run ),
                pybind11::arg( "engine" ), pybind11::arg( "steps" ),
                pybind11::return_value_policy::reference_internal,
                pybind11::call_guard< pybind11::gil_scoped_release >(),
                "Run the given number of Metropolis-Hastings steps from the "
                "current state and return the object sets, which the runner "
                "keeps (they change at the next run)." )

            .def( "run",
                pybind11::overload_cast< geode::RandomEngine&,
                    const geode::SimulationConfigurator& >( &Runner::run ),
                pybind11::arg( "engine" ), pybind11::arg( "config" ),
                // the returned tracker refers to the runner model
                pybind11::keep_alive< 0, 1 >(),
                pybind11::call_guard< pybind11::gil_scoped_release >(),
                "Run the simulation (burn-in, then the realizations) from "
                "the current state, writing the files of config.printer if "
                "any, and return the StatisticsTracker of the realizations." )

            .def( "state_realization", &Runner::state_realization,
                pybind11::return_value_policy::reference_internal,
                "Return the current state of the object sets (the last "
                "realization after a run); it changes at the next run." )

            .def(
                "term_names",
                []( const Runner& self ) {
                    return self.model().term_names();
                },
                "Return the names of the model terms, i.e. of the tracked "
                "statistics." )

            .def(
                "validate_statistics",
                []( const Runner& self,
                    const geode::StatisticsTracker< object_type >& tracker,
                    const std::vector< geode::TargetStatisticConfig >&
                        targets ) {
                    const geode::TargetStatistics< object_type > target_stats{
                        self.model(), targets
                    };
                    geode::statistics::validate( tracker, target_stats );
                },
                pybind11::arg( "tracker" ), pybind11::arg( "targets" ),
                "Check the statistics tracked by a run against the targets "
                "(TargetStatisticConfig, by term name). Raise an "
                "OpenGeodeResultError if the mean of a target statistic has "
                "a relative error above its tolerance, an OpenGeodeDataError "
                "if a term is unknown (or the tracker comes from another "
                "runner)." );
    }

    void define_simulation_configurator_and_runner( pybind11::module_& module )
    {
        define_simulation_configurator( module );
        define_simulation_runner< OwnerSegment2D >( module, "Segment2D" );
    }

} // namespace geode
