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
#include <geode/stochastic/models/energy_term_collection.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
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

        StatisticsTracker< ObjectType > run(
            RandomEngine& engine, const SimulationConfigurator& config )
        {
            if( config.burn_in_steps > 0 )
            {
                mh_sampler_->walk( object_sets_, engine, config.burn_in_steps );
            }

            // Initialize monitoring
            StatisticsTracker< ObjectType > stats_monitor( *model_ );
            std::unique_ptr< SimulationPrinter< ObjectType > > printer;

            if( config.printer.has_value() )
            {
                printer = std::make_unique< SimulationPrinter< ObjectType > >(
                    *model_, config.printer.value() );
            }

            for( const auto realization : Range{ config.realizations } )
            {
                mh_sampler_->walk(
                    object_sets_, engine, config.metropolis_hasting_steps );

                const auto stats = model_->compute_statistics( object_sets_ );
                stats_monitor.add_realization( stats );

                if( printer )
                {
                    printer->print_statistics( stats );
                    printer->print_object_sets( object_sets_, realization );
                }
            }

            if( printer )
            {
                printer->print_statistics_summary( stats_monitor );
            }

            return stats_monitor;
        }

        const ObjectSets< ObjectType >& state_realization() const
        {
            return object_sets_;
        }

    protected:
        // void initialize_sets_and_samplers() = 0;
        // void initialize_model() = 0;

    protected:
        SpatialDomain< ObjectType::dim > domain_;

        ObjectSets< ObjectType > object_sets_;
        std::vector< std::unique_ptr< geode::ObjectSetSampler< ObjectType > > >
            set_samplers_;
        std::unique_ptr< Model< ObjectType > > model_;
        std::unique_ptr< geode::MetropolisHastings< ObjectType > > mh_sampler_;
    };
} // namespace geode