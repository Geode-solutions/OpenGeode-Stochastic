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

#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/stochastic/inference/target_statistics.hpp>

#include <geode/stochastic/models/model.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>

#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < typename ObjectType >
    struct SimulationContext
    {
        std::string string() const
        {
            auto message = std::string{ "SimulationContext: " };
            absl::StrAppend( &message, "\n\t --> ", domain->string() );
            absl::StrAppend( &message, "\n\t --> ", object_sets->string() );
            absl::StrAppend(
                &message, "\n\t --> ", set_samplers.size(), " Sets samplers " );
            absl::StrAppend( &message, "\n\t --> ", model->string() );
            // absl::StrAppend( &message, "\n\t --> ", mh_sampler_ > string() );

            return message;
        }

        std::unique_ptr< SpatialDomain< ObjectType::dim > > domain;

        std::unique_ptr< ObjectSets< ObjectType > > object_sets{
            std::make_unique< ObjectSets< ObjectType > >()
        };
        std::vector< std::unique_ptr< geode::ObjectSetSampler< ObjectType > > >
            set_samplers;
        std::unique_ptr< Model< ObjectType > > model;
        std::unique_ptr< geode::MetropolisHastings< ObjectType > > mh_sampler;

        std::optional< TargetStatistics< ObjectType > > target_statistics;
    };

} // namespace geode