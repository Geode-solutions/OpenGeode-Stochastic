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

#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace geode
{
    void define_simulation_runner( pybind11::module& module )
    {
        pybind11::class_< SimulationConfigurator >(
            module, "SimulationConfigurator" )
            .def( pybind11::init<>() )
            .def_readwrite( "realizations",
                &SimulationConfigurator::realizations,
                "Number of realizations to generate" )
            .def_readwrite( "metropolis_hasting_steps",
                &SimulationConfigurator::metropolis_hasting_steps,
                "Number of Metropolis-Hastings steps per realization" )
            .def_readwrite( "burn_in_steps",
                &SimulationConfigurator::burn_in_steps,
                "Number of burn-in steps before recording realizations" )
            .def_readwrite( "printer", &SimulationConfigurator::printer,
                "Optional SimulationPrinter for output" )
            .def( "string", &SimulationConfigurator::string,
                "Return a detailed description of the object set simulation "
                "configurator" )
            .def( "__repr__", []( const SimulationConfigurator& self ) {
                return "<SimulationConfigurator realizations="
                       + std::to_string( self.realizations ) + ">";
            } );
    }
} // namespace geode