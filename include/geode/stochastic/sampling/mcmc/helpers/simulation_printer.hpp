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

#include <geode/stochastic/common.hpp>

#include <geode/stochastic/inference/statistics_tracker.hpp>
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

    template < typename ObjectType >
    class SimulationPrinter
    {
    public:
        SimulationPrinter( const Model< ObjectType >& model,
            const SimulationPrinterConfigurator& config )
            : model_( model ), config_( config )
        {
        }

        // Print statistics to the configured statistics file
        void print_statistics( const std::vector< double >& stats ) const
        {
            if( !config_.print_statistics )
            {
                return;
            }
            if( !stats_file_path_ )
            {
                stats_file_path_ =
                    ( std::filesystem::path( config_.output_folder )
                        / config_.statistics_filename )
                        .string();
                create_file_and_write_header( *stats_file_path_,
                    absl::StrCat( "# Simulation Statistics\n",
                        energy_terms_name_header() ) );
            }
            std::ofstream file =
                open_file_with_dirs( *stats_file_path_, std::ios::app );
            file << absl::StrJoin( stats, " ; " ) << "\n";
        }

        void print_statistics_summary(
            const StatisticsTracker< ObjectType >& tracker ) const
        {
            if( !config_.print_statistics_summary )
            {
                return;
            }

            const auto summary_path =
                ( std::filesystem::path( config_.output_folder )
                    / config_.statistics_summary_filename )
                    .string();
            create_file_and_write_header(
                summary_path, absl::StrCat( "# Summary statistics\n",
                                  energy_terms_name_header() ) );
            std::ofstream file =
                open_file_with_dirs( summary_path, std::ios::app );
            file << absl::StrCat(
                "# Count:\n", tracker.statiscal_count(), "\n" );
            file << absl::StrCat(
                "# Means:\n", absl::StrJoin( tracker.means(), " ; " ), "\n" );
            file << absl::StrCat( "# Variances:\n",
                absl::StrJoin( tracker.variances(), " ; " ), "\n" );
        }

        void print_object_sets( const ObjectSets< ObjectType >& object_sets,
            index_t realization_id ) const
        {
            if( !config_.print_realisations
                || realization_id % config_.realisations_print_frequency != 0 )
            {
                return;
            }

            // const auto filename =
            //     ( std::filesystem::path( config_.output_folder )
            //         / absl::StrCat(
            //             config_.realisations_prefix, realization_id, ".txt" )
            //             )
            //         .string();
            const auto filename = std::filesystem::path( config_.output_folder )
                                  / absl::StrCat( config_.realisations_prefix,
                                      realization_id, ".txt" );

            std::ofstream file =
                open_file_with_dirs( filename, std::ofstream::out );

            const auto all_objects = object_sets.get_all_object();
            file << "#nb_objects\t" << all_objects.size() << "\n";

            for( const auto& object_id : all_objects )
            {
                const auto& object = object_sets.get_object( object_id );
                file << object.string() << "\t" << object_id.set_id.string()
                     << "\n";
            }
        }

    private:
        void create_file_and_write_header(
            const std::filesystem::path& filename,
            absl::string_view header ) const
        {
            std::ofstream file =
                open_file_with_dirs( filename, std::ofstream::out );

            file << header;
        }

        std::ofstream open_file_with_dirs(
            const std::filesystem::path& file_path,
            std::ios::openmode mode ) const
        {
            auto absolute_path = file_path;
            if( !absolute_path.has_parent_path() )
            {
                absolute_path = std::filesystem::current_path() / absolute_path;
            }
            if( absolute_path.has_parent_path() )
            {
                std::filesystem::create_directories(
                    absolute_path.parent_path() );
            }
            std::ofstream file( absolute_path, mode );
            if( !file.is_open() )
            {
                throw geode::OpenGeodeException(
                    "Cannot open file: " + absolute_path.string() );
            }
            return file;
        }

        std::string energy_terms_name_header() const
        {
            return absl::StrCat(
                absl::StrJoin( model_.term_names(), " ; " ), "\n" );
        }

    private:
        const Model< ObjectType >& model_;
        const SimulationPrinterConfigurator config_;
        mutable std::optional< std::string > stats_file_path_;
    };

} // namespace geode