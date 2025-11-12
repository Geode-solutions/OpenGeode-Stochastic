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

#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>

namespace geode
{
    void define_simulation_printer( pybind11::module& module )
    {
        pybind11::class_< SimulationPrinterConfigurator >(
            module, "SimulationPrinterConfigurator" )
            .def( pybind11::init<>() )
            .def_readwrite( "print_statistics",
                &SimulationPrinterConfigurator::print_statistics )
            .def_readwrite( "statistics_filename",
                &SimulationPrinterConfigurator::statistics_filename )
            .def_readwrite( "print_statistics_summary",
                &SimulationPrinterConfigurator::print_statistics_summary )
            .def_readwrite( "statistics_summary_filename",
                &SimulationPrinterConfigurator::statistics_summary_filename )
            .def_readwrite( "print_realisations",
                &SimulationPrinterConfigurator::print_realisations )
            .def_readwrite( "realisations_prefix",
                &SimulationPrinterConfigurator::realisations_prefix )
            .def_readwrite( "realisations_print_frequency",
                &SimulationPrinterConfigurator::realisations_print_frequency )
            .def_readwrite(
                "output_folder", &SimulationPrinterConfigurator::output_folder )
            .def( "__repr__", []( const SimulationPrinterConfigurator& self ) {
                return "<SimulationPrinterConfigurator output_folder='"
                       + self.output_folder + "'>";
            } );

        //        pybind11::class_< SimulationPrinter >( module,
        //        "SimulationPrinter" )
        //            .def( pybind11::init< const SimulationPrinterConfigurator&
        //            >(),
        //                pybind11::arg( "config" ) )
        //            .def( "print_statistics",
        //            &SimulationPrinter::print_statistics,
        //                pybind11::arg( "stats" ), pybind11::arg( "header" ) =
        //                "", "Print statistics vector to file." )
        //            .def( "print_statistics_summary",
        //                &SimulationPrinter::print_statistics_summary,
        //                pybind11::arg( "monitor" ),
        //                pybind11::arg( "energy_term_names" ) = "",
        //                "Print statistics summary from a StatisticsMonitor."
        //                );
    }
} // namespace geode