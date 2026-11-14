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
    void define_simulation_printer( pybind11::module& module )
    {
        using Config = SimulationPrinterConfigurator;
        pybind11::class_< Config > config( module,
            "SimulationPrinterConfigurator",
            "Output files of a simulation: the statistics of each "
            "realization, their summary and the object sets of some "
            "realizations. The folders are created when needed." );
        config.def(
            pybind11::init(
                []( std::optional< std::string > output_folder,
                    std::optional< bool > print_statistics,
                    std::optional< std::string > statistics_filename,
                    std::optional< bool > print_statistics_summary,
                    std::optional< std::string > statistics_summary_filename,
                    std::optional< bool > print_realisations,
                    std::optional< std::string > realisations_prefix,
                    std::optional< index_t > realisations_print_frequency ) {
                    Config config;
                    detail::set_if_given( config.output_folder, output_folder );
                    detail::set_if_given(
                        config.print_statistics, print_statistics );
                    detail::set_if_given(
                        config.statistics_filename, statistics_filename );
                    detail::set_if_given( config.print_statistics_summary,
                        print_statistics_summary );
                    detail::set_if_given( config.statistics_summary_filename,
                        statistics_summary_filename );
                    detail::set_if_given(
                        config.print_realisations, print_realisations );
                    detail::set_if_given(
                        config.realisations_prefix, realisations_prefix );
                    detail::set_if_given( config.realisations_print_frequency,
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
            "Create a configuration. Defaults: output_folder is the "
            "current directory when the configuration is created, "
            "everything is printed, in statistics.txt, "
            "statistics_summary.txt and pattern_<realization>.txt, one "
            "realization out of 100." );
        def_readwrite_named( config, "print_statistics",
            &Config::print_statistics,
            "Print the statistics of each realization." );
        def_readwrite_named( config, "statistics_filename",
            &Config::statistics_filename,
            "File of the statistics of the realizations, in output_folder." );
        def_readwrite_named( config, "print_statistics_summary",
            &Config::print_statistics_summary,
            "Print the mean and variance of the statistics at the end of the "
            "simulation." );
        def_readwrite_named( config, "statistics_summary_filename",
            &Config::statistics_summary_filename,
            "File of the summary of the statistics, in output_folder." );
        def_readwrite_named( config, "print_realisations",
            &Config::print_realisations,
            "Print the object sets of some realizations." );
        def_readwrite_named( config, "realisations_prefix",
            &Config::realisations_prefix,
            "Prefix of the files of the printed realizations, followed by "
            "the index of the realization and .txt." );
        def_readwrite_named( config, "realisations_print_frequency",
            &Config::realisations_print_frequency,
            "Print one realization out of this number (the realizations 0, "
            "f, 2f...); must be positive when print_realisations is True." );
        def_readwrite_named( config, "output_folder", &Config::output_folder,
            "Folder of the output files." );
        config
            .def( "string", &Config::string,
                "Return a readable description of the configuration." )
            .def( "__repr__", []( const Config& self ) {
                return python_call_repr( "SimulationPrinterConfigurator",
                    { "output_folder=" + python_repr( self.output_folder ),
                        "print_statistics="
                            + python_repr( self.print_statistics ),
                        "statistics_filename="
                            + python_repr( self.statistics_filename ),
                        "print_statistics_summary="
                            + python_repr( self.print_statistics_summary ),
                        "statistics_summary_filename="
                            + python_repr( self.statistics_summary_filename ),
                        "print_realisations="
                            + python_repr( self.print_realisations ),
                        "realisations_prefix="
                            + python_repr( self.realisations_prefix ),
                        "realisations_print_frequency="
                            + python_repr(
                                self.realisations_print_frequency ) } );
            } );
    }
} // namespace geode
