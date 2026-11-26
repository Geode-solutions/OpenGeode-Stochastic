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
#include <filesystem>
#include <fstream>
#include <string>

#include <absl/strings/match.h>
#include <absl/strings/str_join.h>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/applications/strauss_process.hpp>
#include <geode/stochastic/inference/statistics_tracker.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>

namespace
{
    // Model with two terms: the density of the set and the interaction
    geode::SimulationContext< geode::Point2D > two_term_context()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { 10., 10. } };
        domain.buffer_size = 1.;
        geode::StraussProcessBuilder< geode::Point2D > strauss;
        strauss.set_domain( domain );
        auto& set = strauss.add_set( "points", 0.5 );
        geode_unused( set );
        strauss.add_interaction( { "points" }, 0.5, 1., std::nullopt );
        // NOLINTEND(*-magic-numbers)
        return strauss.build_simulation_context();
    }

    std::vector< std::string > read_lines( const std::filesystem::path& path )
    {
        geode::OpenGeodeStochasticStochasticException::test(
            std::filesystem::exists( path ),
            "[SimulationPrinter] File not created: ", path.string() );
        std::ifstream file{ path };
        std::vector< std::string > lines;
        std::string line;
        while( std::getline( file, line ) )
        {
            lines.push_back( line );
        }
        return lines;
    }

    void test_print_statistics( const geode::Model< geode::Point2D >& model,
        const geode::SimulationPrinterConfigurator& config )
    {
        geode::Logger::info( "TEST - SIMULATION PRINTER STATISTICS" );
        const geode::SimulationPrinter< geode::Point2D > printer{ model,
            config };
        // NOLINTBEGIN(*-magic-numbers)
        printer.print_statistics( { 1., 2.5 } );
        printer.print_statistics( { 3., 4. } );
        // NOLINTEND(*-magic-numbers)

        const auto lines =
            read_lines( std::filesystem::path{ config.output_folder }
                        / config.statistics_filename );
        const std::vector< std::string > expected{ "# Simulation Statistics",
            absl::StrJoin( model.term_names(), " ; " ), "1 ; 2.5", "3 ; 4" };
        geode::OpenGeodeStochasticStochasticException::test( lines == expected,
            "[SimulationPrinter] Wrong statistics file: ",
            absl::StrJoin( lines, "\n" ) );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_print_statistics_summary(
        const geode::Model< geode::Point2D >& model,
        const geode::SimulationPrinterConfigurator& config )
    {
        geode::Logger::info( "TEST - SIMULATION PRINTER SUMMARY" );
        const geode::SimulationPrinter< geode::Point2D > printer{ model,
            config };
        geode::StatisticsTracker< geode::Point2D > tracker{ model };
        // NOLINTBEGIN(*-magic-numbers)
        tracker.add_realization( { 1., 2. } );
        tracker.add_realization( { 3., 2. } );
        // NOLINTEND(*-magic-numbers)
        printer.print_statistics_summary( tracker );

        const auto lines =
            read_lines( std::filesystem::path{ config.output_folder }
                        / config.statistics_summary_filename );
        const std::vector< std::string > expected{ "# Summary statistics",
            absl::StrJoin( model.term_names(), " ; " ), "# Count:", "2",
            "# Means:", "2 ; 2", "# Variances:", "2 ; 0" };
        geode::OpenGeodeStochasticStochasticException::test( lines == expected,
            "[SimulationPrinter] Wrong summary file: ",
            absl::StrJoin( lines, "\n" ) );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_print_object_sets( const geode::Model< geode::Point2D >& model,
        geode::SimulationPrinterConfigurator config )
    {
        geode::Logger::info( "TEST - SIMULATION PRINTER OBJECT SETS" );
        geode::ObjectSets< geode::Point2D > object_sets;
        const auto set_id = object_sets.add_set( "points" );
        // NOLINTBEGIN(*-magic-numbers)
        object_sets.add_object( geode::Point2D{ { 0., 0. } }, set_id, false );
        object_sets.add_object( geode::Point2D{ { 1., 0. } }, set_id, false );
        object_sets.add_object( geode::Point2D{ { 3., 0. } }, set_id, true );
        config.realisations_print_frequency = 2;
        // NOLINTEND(*-magic-numbers)
        const geode::SimulationPrinter< geode::Point2D > printer{ model,
            config };
        for( const auto realization : geode::Range{ 3 } )
        {
            printer.print_object_sets( object_sets, realization );
        }

        const std::filesystem::path folder{ config.output_folder };
        const auto prefix = config.realisations_prefix;
        geode::OpenGeodeStochasticStochasticException::test(
            !std::filesystem::exists( folder / ( prefix + "1.txt" ) ),
            "[SimulationPrinter] A realization out of the print frequency "
            "should not be printed." );
        geode::OpenGeodeStochasticStochasticException::test(
            std::filesystem::exists( folder / ( prefix + "2.txt" ) ),
            "[SimulationPrinter] Realization 2 should be printed." );
        const auto lines = read_lines( folder / ( prefix + "0.txt" ) );
        geode::OpenGeodeStochasticStochasticException::test(
            lines.size() == 4 && lines[0] == "#nb_objects\t3",
            "[SimulationPrinter] Wrong realization file: ",
            absl::StrJoin( lines, "\n" ) );
        for( const auto line : geode::Range{ 1, 4 } )
        {
            geode::OpenGeodeStochasticStochasticException::test(
                absl::EndsWith( lines[line], "\tpoints" ),
                "[SimulationPrinter] Each object line should end with its set "
                "name: ",
                lines[line] );
        }

        geode::Logger::info( "--> SUCCESS!" );
    }

    // A new printer restarts the statistics file, and a summary replaces the
    // previous one: SimulationRunner::run starts a new series, with a new
    // printer, so two runs in the same folder do not mix their outputs.
    // Extending a series is a separate feature (T17 of the analysis).
    void test_files_are_rewritten( const geode::Model< geode::Point2D >& model,
        const geode::SimulationPrinterConfigurator& config )
    {
        geode::Logger::info( "TEST - SIMULATION PRINTER REWRITES ITS FILES" );
        const std::filesystem::path folder{ config.output_folder };
        // NOLINTBEGIN(*-magic-numbers)
        const geode::SimulationPrinter< geode::Point2D > printer{ model,
            config };
        printer.print_statistics( { 5., 6. } );
        const auto stats_lines =
            read_lines( folder / config.statistics_filename );
        geode::OpenGeodeStochasticStochasticException::test(
            stats_lines.size() == 3 && stats_lines[2] == "5 ; 6",
            "[SimulationPrinter] A new printer should restart the statistics "
            "file: ",
            absl::StrJoin( stats_lines, "\n" ) );

        geode::StatisticsTracker< geode::Point2D > tracker{ model };
        tracker.add_realization( { 1., 2. } );
        printer.print_statistics_summary( tracker );
        tracker.add_realization( { 3., 4. } );
        printer.print_statistics_summary( tracker );
        // NOLINTEND(*-magic-numbers)
        const auto summary_lines =
            read_lines( folder / config.statistics_summary_filename );
        geode::OpenGeodeStochasticStochasticException::test(
            summary_lines.size() == 8 && summary_lines[3] == "2",
            "[SimulationPrinter] The summary should replace the previous one: ",
            absl::StrJoin( summary_lines, "\n" ) );

        geode::Logger::info( "--> SUCCESS!" );
    }

    // Nothing is written when the outputs are off; the print frequency is
    // checked only when the realizations are printed.
    void test_disabled_outputs( const geode::Model< geode::Point2D >& model,
        geode::SimulationPrinterConfigurator config )
    {
        geode::Logger::info( "TEST - SIMULATION PRINTER DISABLED OUTPUTS" );
        config.output_folder =
            ( std::filesystem::path{ config.output_folder } / "disabled" )
                .string();
        config.realisations_print_frequency = 0;
        bool rejected{ false };
        try
        {
            const geode::SimulationPrinter< geode::Point2D > printer{ model,
                config };
        }
        catch( const geode::OpenGeodeException& )
        {
            rejected = true;
        }
        geode::OpenGeodeStochasticStochasticException::test( rejected,
            "[SimulationPrinter] A zero print frequency should be rejected." );

        config.print_statistics = false;
        config.print_statistics_summary = false;
        config.print_realisations = false;
        const geode::SimulationPrinter< geode::Point2D > printer{ model,
            config };
        geode::StatisticsTracker< geode::Point2D > tracker{ model };
        // NOLINTBEGIN(*-magic-numbers)
        tracker.add_realization( { 1., 2. } );
        printer.print_statistics( { 1., 2. } );
        // NOLINTEND(*-magic-numbers)
        printer.print_statistics_summary( tracker );
        geode::ObjectSets< geode::Point2D > object_sets;
        object_sets.add_set( "points" );
        printer.print_object_sets( object_sets, 0 );
        geode::OpenGeodeStochasticStochasticException::test(
            !std::filesystem::exists( config.output_folder ),
            "[SimulationPrinter] Nothing should be written when the outputs "
            "are off." );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        const auto context = two_term_context();
        const auto folder =
            std::filesystem::current_path() / "simulation_printer_test";
        std::filesystem::remove_all( folder );

        geode::SimulationPrinterConfigurator config;
        config.output_folder = folder.string();
        config.statistics_filename = "stats.txt";
        config.statistics_summary_filename = "summary.txt";
        config.realisations_prefix = "pattern_";

        test_print_statistics( *context.model, config );
        test_print_statistics_summary( *context.model, config );
        test_print_object_sets( *context.model, config );
        test_files_are_rewritten( *context.model, config );
        test_disabled_outputs( *context.model, config );
        std::filesystem::remove_all( folder );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
