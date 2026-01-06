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

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_monitor.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>

#include <filesystem>
#include <fstream>

namespace
{
    void test_print_statistics(
        const geode::SimulationPrinterConfigurator &config )
    {
        geode::Logger::info( "[TEST] SimulationPrinter print statistics" );

        geode::SimulationPrinter printer( config );

        // --- Test print_statistics
        std::vector< double > stats = { 1.0, 2.5, 3.7 };
        printer.print_statistics(
            stats, "EnergyTerm1;EnergyTerm2;EnergyTerm3" );

        const std::filesystem::path temp_folder = config.output_folder;
        const auto stats_path = temp_folder / config.statistics_filename;
        OPENGEODE_EXCEPTION( std::filesystem::exists( stats_path ),
            "Statistics file not created" );

        std::ifstream stats_file( stats_path );
        std::string line;
        std::getline( stats_file, line );
        OPENGEODE_EXCEPTION(
            line == "# Simulation Statistics", "Header not correctly written" );
        std::getline( stats_file, line );
        OPENGEODE_EXCEPTION( line == "EnergyTerm1;EnergyTerm2;EnergyTerm3",
            "Header not correctly written" );
        std::getline( stats_file, line );
        OPENGEODE_EXCEPTION(
            line == "1 ; 2.5 ; 3.7", "Statistics line not correctly written" );

        geode::Logger::info( "--> Success!" );
    }

    void test_statistics_summary(
        const geode::SimulationPrinterConfigurator &config )
    {
        geode::Logger::info(
            "[TEST] SimulationPrinter print statistics summary" );

        geode::SimulationPrinter printer( config );
        geode::StatisticsMonitor monitor( 2 );
        monitor.add_realization( { 1, 2 } );

        printer.print_statistics_summary( monitor, "EnergyTerm1;EnergyTerm2" );

        const std::filesystem::path temp_folder = config.output_folder;
        const auto summary_path =
            temp_folder / config.statistics_summary_filename;

        OPENGEODE_EXCEPTION( std::filesystem::exists( summary_path ),
            "Summary file not created" );

        std::ifstream summary_file( summary_path );
        std::string content(
            ( std::istreambuf_iterator< char >( summary_file ) ),
            std::istreambuf_iterator< char >() );
        OPENGEODE_EXCEPTION(
            content.find( "EnergyTerm1;EnergyTerm2" ) != std::string::npos,
            "Energy term names missing" );
        OPENGEODE_EXCEPTION(
            content.find( "2" ) != std::string::npos, "Count missing" );
        OPENGEODE_EXCEPTION(
            content.find( "1 ; 2" ) != std::string::npos, "Means missing" );
        OPENGEODE_EXCEPTION(
            content.find( "0 ; 0" ) != std::string::npos, "Variances missing" );

        geode::Logger::info( "--> Success!" );
    }

    void test_print_objects(
        const geode::SimulationPrinterConfigurator &config )
    {
        geode::Logger::info( "[TEST] SimulationPrinter print object" );

        geode::ObjectSets< geode::Point2D > object_sets;
        const auto set_id = object_sets.add_set( "default_name" );

        object_sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id, false );
        object_sets.add_object( geode::Point2D{ { 1.0, 0.0 } }, set_id, false );
        object_sets.add_object( geode::Point2D{ { 3.0, 0.0 } }, set_id, false );

        geode::SimulationPrinter printer( config );
        printer.print_object_sets( object_sets, 0 );

        const std::filesystem::path temp_folder = config.output_folder;

        const auto obj_path = temp_folder / "pattern_0.txt";
        OPENGEODE_EXCEPTION( std::filesystem::exists( obj_path ),
            "Object sets file not created" );

        geode::Logger::info( "--> Success!" );
    }

} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        const std::filesystem::path temp_folder =
            std::filesystem::current_path() / "simprinter_test";
        SDEBUG( temp_folder );
        std::filesystem::create_directories( temp_folder );

        geode::SimulationPrinterConfigurator config;
        config.output_folder = temp_folder.string();
        config.statistics_filename = "stats.txt";
        config.statistics_summary_filename = "summary.txt";
        config.realisations_prefix = "pattern_";

        test_print_statistics( config );
        test_statistics_summary( config );
        test_print_objects( config );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}