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

#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>

namespace geode
{
    void define_simulation_printer( pybind11::module &module )
    {
        pybind11::class_< SimulationPrinterConfigurator >(
            module, "SimulationPrinterConfigurator" )
            .def( pybind11::init(
                      []( std::optional< std::string > output_folder,
                          std::optional< bool > print_statistics,
                          std::optional< std::string > statistics_filename,
                          std::optional< bool > print_statistics_summary,
                          std::optional< std::string >
                              statistics_summary_filename,
                          std::optional< bool > print_realisations,
                          std::optional< std::string > realisations_prefix,
                          std::optional< index_t >
                              realisations_print_frequency ) {
                          SimulationPrinterConfigurator config;
                          detail::set_if_given(
                              config.output_folder, output_folder );
                          detail::set_if_given(
                              config.print_statistics, print_statistics );
                          detail::set_if_given(
                              config.statistics_filename, statistics_filename );
                          detail::set_if_given( config.print_statistics_summary,
                              print_statistics_summary );
                          detail::set_if_given(
                              config.statistics_summary_filename,
                              statistics_summary_filename );
                          detail::set_if_given(
                              config.print_realisations, print_realisations );
                          detail::set_if_given(
                              config.realisations_prefix, realisations_prefix );
                          detail::set_if_given(
                              config.realisations_print_frequency,
                              realisations_print_frequency );
                          return config;
                      } ),
                pybind11::kw_only(),
                pybind11::arg( "output_folder" ) = std::nullopt,
                pybind11::arg( "print_statistics" ) = std::nullopt,
                pybind11::arg( "statistics_filename" ) = std::nullopt,
                pybind11::arg( "print_statistics_summary" ) = std::nullopt,
                pybind11::arg( "statistics_summary_filename" ) = std::nullopt,
                pybind11::arg( "print_realisations" ) = std::nullopt,
                pybind11::arg( "realisations_prefix" ) = std::nullopt,
                pybind11::arg( "realisations_print_frequency" ) = std::nullopt,
                "Output files of a simulation; output_folder defaults to the "
                "current directory." )
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
            .def( "__repr__", []( const SimulationPrinterConfigurator &self ) {
                return "<SimulationPrinterConfigurator output_folder='"
                       + self.output_folder + "'>";
            } );
    }
} // namespace geode