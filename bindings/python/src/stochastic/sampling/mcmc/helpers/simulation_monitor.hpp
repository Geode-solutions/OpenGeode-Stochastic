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

#include <geode/stochastic/sampling/mcmc/helpers/simulation_monitor.hpp>

namespace geode
{
    void define_simulation_monitor( pybind11::module &module )
    {
        pybind11::class_< geode::StatisticsMonitor >(
            module, "StatisticsMonitor" )
            .def( pybind11::init< geode::index_t >(),
                pybind11::arg( "nb_energy_terms" ),
                "Create a StatisticsMonitor for a given number of energy "
                "terms" )
            .def( "add_realization", &geode::StatisticsMonitor::add_realization,
                pybind11::arg( "values" ),
                "Add a realization (vector of doubles) to update statistics" )
            .def( "statiscal_count", &geode::StatisticsMonitor::statiscal_count,
                "Return the number of realizations added" )
            .def_property_readonly( "means", &geode::StatisticsMonitor::means,
                "Return the computed mean values for each energy term" )
            .def_property_readonly( "variances",
                &geode::StatisticsMonitor::variances,
                "Return the computed variances for each energy term" )
            .def( "__repr__", []( const geode::StatisticsMonitor &self ) {
                return "<StatisticsMonitor count="
                       + std::to_string( self.statiscal_count() ) + ">";
            } );
    }
} // namespace geode