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
#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/energy_term_collection.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

#include <absl/strings/str_join.h>

namespace geode
{
    struct SimulationConfigurator
    {
        index_t realizations{ 1000 };
        index_t metropolis_hasting_steps{ 1000 };
        index_t burn_in_steps{ 1000 };

        std::optional< SimulationPrinterConfigurator > printer{ std::nullopt };

        std::string string() const
        {
            auto message = absl::StrCat( "SimulationConfigurator: " );
            absl::StrAppend( &message, "\n\t --> ", realizations,
                " metropolis hasting realizations" );
            absl::StrAppend( &message, "\n\t --> ", metropolis_hasting_steps,
                " metropolis hasting steps" );
            absl::StrAppend(
                &message, "\n\t --> ", burn_in_steps, " burnin steps" );
            if( printer.has_value() )
            {
                absl::StrAppend( &message, "\n\t --> Simulation Printer: \n",
                    printer.value().string() );
            }
            return message;
        }
    };

    template < typename ObjectType >
    class SimulationRunner
    {
    public:
        SimulationRunner( const SpatialDomain< ObjectType::dim >& domain )
            : domain_( domain ) {};
        virtual ~SimulationRunner() = default;

        virtual void initialize() = 0;

        const ObjectSets< ObjectType >& run(
            RandomEngine& engine, const index_t steps )
        {
            mh_sampler_->walk( object_sets_, engine, steps );
            return object_sets_;
        }

        StatisticsMonitor run(
            RandomEngine& engine, const SimulationConfigurator& config )
        {
            if( config.burn_in_steps > 0 )
            {
                mh_sampler_->walk( object_sets_, engine, config.burn_in_steps );
            }

            // Initialize monitoring
            StatisticsMonitor stats_monitor( energy_terms_collection_.size() );
            std::unique_ptr< SimulationPrinter > printer;

            if( config.printer.has_value() )
            {
                printer = std::make_unique< SimulationPrinter >(
                    config.printer.value() );
            }

            for( const auto realization : Range{ config.realizations } )
            {
                mh_sampler_->walk(
                    object_sets_, engine, config.metropolis_hasting_steps );

                const auto stats = state_statistics();
                stats_monitor.add_realization( stats );

                if( printer )
                {
                    printer->print_statistics(
                        stats, model_energy_term_names() );
                    printer->print_object_sets( object_sets_, realization );
                }
            }

            if( printer )
            {
                printer->print_statistics_summary(
                    stats_monitor, model_energy_term_names() );
            }

            return stats_monitor;
        }

        const ObjectSets< ObjectType >& state_realization() const
        {
            return object_sets_;
        }

        std::vector< double > state_statistics() const
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

        std::string model_energy_term_names() const
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

    protected:
        SpatialDomain< ObjectType::dim > domain_;
        std::vector< std::unique_ptr< geode::ObjectSetSampler< ObjectType > > >
            set_samplers_;

        std::vector< geode::uuid > ordered_energy_terms_;
        std::vector< double > ordered_target_statistics_;

        EnergyTermCollection< ObjectType > energy_terms_collection_;
        std::unique_ptr< geode::MetropolisHastings< ObjectType > > mh_sampler_;

        ObjectSets< ObjectType > object_sets_;
    };
} // namespace geode