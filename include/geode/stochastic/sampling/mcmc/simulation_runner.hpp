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
#include <geode/stochastic/sampling/mcmc/helpers/simulation_printer.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>

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
        SimulationRunner( SimulationContext< ObjectType >&& context )
            : context_( std::move( context ) ) {};
        virtual ~SimulationRunner() = default;

        const ObjectSets< ObjectType >& run(
            RandomEngine& engine, const index_t steps )
        {
            context_.mh_sampler->walk( *context_.object_sets, engine, steps );
            return *context_.object_sets;
        }

        StatisticsTracker< ObjectType > run(
            RandomEngine& engine, const SimulationConfigurator& config )
        {
            if( config.burn_in_steps > 0 )
            {
                context_.mh_sampler->walk(
                    *context_.object_sets, engine, config.burn_in_steps );
            }

            // Initialize monitoring
            StatisticsTracker< ObjectType > stats_monitor( *context_.model );
            std::unique_ptr< SimulationPrinter< ObjectType > > printer;

            if( config.printer.has_value() )
            {
                printer = std::make_unique< SimulationPrinter< ObjectType > >(
                    *context_.model, config.printer.value() );
            }

            for( const auto realization : Range{ config.realizations } )
            {
                context_.mh_sampler->walk( *context_.object_sets, engine,
                    config.metropolis_hasting_steps );

                const auto stats =
                    context_.model->compute_statistics( *context_.object_sets );
                stats_monitor.add_realization( stats );

                if( printer )
                {
                    printer->print_statistics( stats );
                    printer->print_object_sets(
                        *context_.object_sets, realization );
                }
            }

            if( printer )
            {
                printer->print_statistics_summary( stats_monitor );
            }

            return stats_monitor;
        }

        [[nodiscard]] const TargetStatistics< ObjectType >&
            target_statistics() const
        {
            OpenGeodeStochasticStochasticException::check_exception(
                context_.target_statistics.has_value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[SimulationRunner] Target statistics not initialized" );

            return *context_.target_statistics;
        }

        [[nodiscard]] const ObjectSets< ObjectType >& state_realization() const
        {
            return *context_.object_sets;
        }

    protected:
        SimulationContext< ObjectType > context_;
    };
} // namespace geode