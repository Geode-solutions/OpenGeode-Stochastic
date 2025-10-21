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
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/energy_term_collection.hpp>

#include <absl/strings/str_join.h>
#include <fstream>

namespace geode
{
    class opengeode_stochastic_stochastic_api MonitoringStatistics
    {
    public:
        std::vector< double > sum;
        std::vector< double > sum_squares;
        std::vector< double > means;
        std::vector< double > variances;
        MonitoringStatistics( MonitoringStatistics&& ) = default;
        MonitoringStatistics( const MonitoringStatistics& ) = default;

        MonitoringStatistics( const index_t nb_energy_terms )
        {
            sum.resize( nb_energy_terms, 0.0 );
            sum_squares.resize( nb_energy_terms, 0.0 );
            means.resize( nb_energy_terms, 0.0 );
            variances.resize( nb_energy_terms, 0.0 );
        }

        void add_realization( const std::vector< double >& values )
        {
            for( const auto stat_id : Range{ values.size() } )
            {
                sum[stat_id] += values[stat_id];
                sum_squares[stat_id] += values[stat_id] * values[stat_id];
            }
        }

        void finalize( const index_t nb_realizations )
        {
            for( const auto stat_id : Range{ sum.size() } )
            {
                means[stat_id] = sum[stat_id] / nb_realizations;
                double variance =
                    ( sum_squares[stat_id]
                        - ( sum[stat_id] * sum[stat_id] ) / nb_realizations )
                    / ( nb_realizations - 1 );
                variances[stat_id] = variance;
                // stddevs[stat_id] =std::sqrt( std::max( variance, 0.0 ) );
            }
        }
    };

    template < typename ObjectType >
    class SimulationRunner
    {
    public:
        SimulationRunner() = default;
        virtual ~SimulationRunner() = default;

        virtual void initialize() = 0;

        const ObjectSets< ObjectType >& run(
            RandomEngine& engine, const index_t steps )
        {
            mh_sampler_->walk( object_sets_, engine, steps );
            return object_sets_;
        }

        void run_and_print( std::string_view filename,
            RandomEngine& engine,
            const index_t steps,
            const index_t nb_realizations )
        {
            const auto file_exist =
                static_cast< bool >( std::ifstream( filename.data() ) );
            if( !file_exist )
            {
                const auto header = statistics_header_file();
                print_to_file( filename, header );
            }

            for( const auto realization : Range{ nb_realizations } )
            {
                run( engine, steps );
                const auto statistics = statistics_string();
                print_to_file( filename, statistics );
            }
        }

        MonitoringStatistics run_print_and_monitor( std::string_view filename,
            RandomEngine& engine,
            const index_t steps,
            const index_t nb_realizations )
        {
            const auto file_exist =
                static_cast< bool >( std::ifstream( filename.data() ) );
            if( !file_exist )
            {
                const auto header = statistics_header_file();
                print_to_file( filename, header );
            }
            MonitoringStatistics stat_monitoring(
                energy_terms_collection_.size() );

            for( const auto realization : Range{ nb_realizations } )
            {
                run( engine, steps );
                const auto stats = statistics();
                print_to_file( filename,
                    absl::StrCat( absl::StrJoin( stats, " ; " ), "\n" ) );
                stat_monitoring.add_realization( stats );
            }
            stat_monitoring.finalize( nb_realizations );
            return stat_monitoring;
        }

        const ObjectSets< ObjectType >& current_pattern_realization() const
        {
            return object_sets_;
        }

        std::vector< double > statistics() const
        {
            std::vector< double > statistic_values;
            statistic_values.reserve( ordered_energy_terms_.size() );

            for( const auto& energy_term_uuid : ordered_energy_terms_ )
            {
                const auto& term =
                    energy_terms_collection_.get( energy_term_uuid );
                statistic_values.push_back( term.statistic( object_sets_ ) );
            }

            return statistic_values;
        }

        std::string statistics_log_info() const
        {
            std::string message( "Pattern statistics: " );
            for( const auto term_id :
                geode::Range{ ordered_energy_terms_.size() } )
            {
                const auto& energy_term = energy_terms_collection_.get(
                    ordered_energy_terms_[term_id] );
                const double value = energy_term.statistic( object_sets_ );
                absl::StrAppend( &message, " \t Term(", energy_term.name(),
                    ") --> value/traget: ", value, " / ",
                    ordered_target_statistics_[term_id] );
            }
            return message;
        }

    protected:
        std::string energy_term_names() const
        {
            std::vector< std::string > term_names;
            term_names.reserve( ordered_energy_terms_.size() );

            for( const auto& energy_term_uuid : ordered_energy_terms_ )
            {
                const auto& term =
                    energy_terms_collection_.get( energy_term_uuid );
                term_names.push_back( term.name().data() );
            }

            return absl::StrCat( absl::StrJoin( term_names, " ; " ), "\n" );
        }

        std::string statistics_string() const
        {
            return absl::StrCat( absl::StrJoin( statistics(), " ; " ), "\n" );
        }

        std::string statistics_header_file()
        {
            std::string message( "Sufficient statistics mcmc iterations:\n" );
            absl::StrAppend( &message, energy_term_names() );
            return message;
        }

        void print_to_file(
            absl::string_view filename, absl::string_view message )
        {
            std::ofstream file(
                filename.data(), std::ofstream::out | std::ofstream::app );
            file << message;
            file.close();
            return;
        }

    protected:
        std::vector< std::unique_ptr< geode::ObjectSetSampler< ObjectType > > >
            set_samplers_;

        std::vector< geode::uuid > ordered_energy_terms_;
        std::vector< double > ordered_target_statistics_;

        EnergyTermCollection< ObjectType > energy_terms_collection_;
        std::unique_ptr< geode::MetropolisHastings< ObjectType > > mh_sampler_;

        ObjectSets< ObjectType > object_sets_;
    };
} // namespace geode