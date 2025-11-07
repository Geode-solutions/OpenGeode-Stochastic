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
// #include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
// #include <geode/stochastic/sampling/mcmc/models/energy_term_collection.hpp>

#include <absl/strings/str_join.h>
#include <filesystem>
#include <fstream>

namespace geode
{
    struct SimulationPrinterConfigurator
    {
        bool print_statistics{ true };
        std::string statistics_filename{ "statistics.txt" };

        bool print_realisations{ true };
        std::string realisations_prefix{ "pattern_" };
        index_t realisations_print_frequency{ 100 };

        std::string output_folder{ "." };
    };

    class SimulationPrinter
    {
    public:
        SimulationPrinter( absl::string_view stats_file )
            : stats_file_( stats_file )
        {
            write_header_if_new( stats_file_, "Simulation Statistics\n" );
        }

        // Print header if file does not exist
        void write_header_if_new(
            absl::string_view filename, absl::string_view header )
        {
            namespace fs = std::filesystem;
            fs::path file_path{ std::string( filename ) };
            if( !fs::exists( file_path ) )
            {
                std::ofstream file = open_file_with_dirs( filename );
                file << header;
            }
        }

        void print_statistics( const std::vector< double >& stats )
        {
            std::ofstream file =
                open_file_with_dirs( stats_file_, std::ios::app );
            file << absl::StrJoin( stats, " ; " ) << "\n";
        }

        template < typename ObjectType >
        void print_object_sets( const ObjectSets< ObjectType >& object_sets,
            absl::string_view filename )
        {
            std::ofstream file = open_file_with_dirs( filename );

            const auto all_objects = object_sets.get_all_object();

            file << "#nb_objects\t" << all_objects.size() << "\n";

            for( const auto& object_id : all_objects )
            {
                const auto& object = object_sets.get_object( object_id );
                file << object.string() << "\t" << object_id.set_id.string()
                     << "\n";
            }
        }

        // Reusable helper to open files with directories created
        std::ofstream open_file_with_dirs( absl::string_view path_filename,
            std::ios::openmode mode = std::ofstream::out )
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

    private:
        std::string stats_file_;
    };

} // namespace geode