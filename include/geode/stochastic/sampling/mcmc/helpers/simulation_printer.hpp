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

#pragma once

#include <geode/stochastic/common.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_monitor.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

#include <absl/strings/str_join.h>
#include <filesystem>
#include <fstream>

namespace geode
{
    struct SimulationPrinterConfigurator
    {
        bool print_statistics{ true };
        std::string statistics_filename{ "statistics.txt" };

        bool print_statistics_summary{ true };
        std::string statistics_summary_filename{ "statistics_summary.txt" };

        bool print_realisations{ true };
        std::string realisations_prefix{ "pattern_" };
        index_t realisations_print_frequency{ 100 };

        std::string output_folder{ std::filesystem::current_path().string() };

        std::string string() const
        {
            auto message = absl::StrCat(
                "SimulationPrinterConfigurator - print to folder: ",
                output_folder );
            if( print_statistics )
            {
                absl::StrAppend( &message,
                    "\n\t --> print statistics for each realisation to file: ",
                    statistics_filename );
            }
            if( print_statistics_summary )
            {
                absl::StrAppend( &message,
                    "\n\t --> print summary statistics to file: ",
                    statistics_summary_filename );
            }
            if( print_realisations )
            {
                absl::StrAppend( &message, "\n\t --> print pattern every ",
                    realisations_print_frequency,
                    " realizations behind the prefix: ", realisations_prefix );
            }
            return message;
        }
    };

    class SimulationPrinter
    {
    public:
        SimulationPrinter( const SimulationPrinterConfigurator &config )
            : config_( config )
        {
        }

        // Print statistics to the configured statistics file
        void print_statistics(
            const std::vector< double > &stats, absl::string_view header ) const
        {
            if( !config_.print_statistics )
                return;
            const auto stats_file_path =
                stats_file_path_.value_or( create_statistics_file( header ) );
            std::ofstream file =
                open_file_with_dirs( stats_file_path, std::ios::app );
            file << absl::StrJoin( stats, " ; " ) << "\n";
        }

        template < typename ObjectType >
        void print_object_sets( const ObjectSets< ObjectType > &object_sets,
            index_t realization_id ) const
        {
            if( !config_.print_realisations
                || realization_id % config_.realisations_print_frequency != 0 )
                return;

            const auto filename =
                ( std::filesystem::path( config_.output_folder )
                    / absl::StrCat(
                        config_.realisations_prefix, realization_id, ".txt" ) )
                    .string();

            std::ofstream file = open_file_with_dirs( filename );

            const auto all_objects = object_sets.get_all_object();
            file << "#nb_objects\t" << all_objects.size() << "\n";

            for( const auto &object_id : all_objects )
            {
                const auto &object = object_sets.get_object( object_id );
                file << object.string() << "\t" << object_id.set_id.string()
                     << "\n";
            }
        }
        void print_statistics_summary( const StatisticsMonitor &monitor,
            absl::string_view energy_term_names ) const
        {
            if( !config_.print_statistics_summary )
                return;

            const auto summary_path =
                ( std::filesystem::path( config_.output_folder )
                    / config_.statistics_summary_filename )
                    .string();

            std::ofstream file = open_file_with_dirs( summary_path );
            file << "# Summary statistics\n";
            file << energy_term_names.data() << "\n";
            file << absl::StrCat(
                "# Count:\n", monitor.statiscal_count(), "\n" );
            file << absl::StrCat(
                "# Means:\n", absl::StrJoin( monitor.means(), " ; " ), "\n" );
            file << absl::StrCat( "# Variances:\n",
                absl::StrJoin( monitor.variances(), " ; " ), "\n" );
        }

    private:
        void write_header_if_new(
            absl::string_view filename, absl::string_view header ) const
        {
            namespace fs = std::filesystem;
            fs::path file_path{ std::string( filename ) };
            if( !fs::exists( file_path ) )
            {
                std::ofstream file = open_file_with_dirs( filename );
                file << header;
            }
        }

        std::ofstream open_file_with_dirs( absl::string_view path_filename,
            std::ios::openmode mode = std::ofstream::out ) const
        {
            namespace fs = std::filesystem;
            fs::path file_path{ std::string( path_filename ) };

            if( !file_path.has_parent_path() )
                file_path = fs::current_path() / file_path;

            if( file_path.has_parent_path() )
                fs::create_directories( file_path.parent_path() );

            std::ofstream file( file_path, mode );
            if( !file.is_open() )
                throw geode::OpenGeodeException(
                    "Cannot open file: " + file_path.string() );

            return file;
        }
        const std::string &create_statistics_file(
            absl::string_view header ) const
        {
            stats_file_path_ = ( std::filesystem::path( config_.output_folder )
                                 / config_.statistics_filename )
                                   .string();

            if( config_.print_statistics )
            {
                write_header_if_new( *stats_file_path_,
                    absl::StrCat(
                        "# Simulation Statistics\n", header.data(), "\n" ) );
            }
            return *stats_file_path_;
        }

    private:
        SimulationPrinterConfigurator config_;
        mutable std::optional< std::string > stats_file_path_;
    };

} // namespace geode