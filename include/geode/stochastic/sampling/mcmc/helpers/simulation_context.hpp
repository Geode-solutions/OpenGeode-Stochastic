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
#include <geode/stochastic/spatial/spatial_domain.hpp>

#include <geode/stochastic/inference/target_statistics.hpp>

#include <geode/stochastic/models/model.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>

#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>

#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/object_set_dynamic_config.hpp>

namespace geode
{
    template < typename ObjectType >
    struct SimulationContext
    {
        [[nodiscard]] std::string string() const
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
    };

    template < typename ObjectType >
    struct SimulationContextConfig
    {
        SpatialDomainConfig< ObjectType::dim > domain;

        std::vector< ObjectSetConfig > sets;
        std::vector< ObjectSetDynamicsConfig > proposals;

        geode::ModelConfig model;
    };

    template < typename ObjectType >
    [[nodiscard]] geode::SimulationContext< ObjectType >
        build_simulation_context(
            const SimulationContextConfig< ObjectType >& config )
    {
        geode::SimulationContext< ObjectType > context;

        // -------------------------
        // Domain
        // -------------------------
        context.domain = geode::build_spatial_domain( config.domain );

        // -------------------------
        // Sets
        // -------------------------

        //        auto proposal_kernel =
        //            std::make_unique< geode::ProposalKernel< geode::Point2D >
        //            >();
        //        for( const auto& set_desc : set_descriptors_ )
        //        {
        //            const auto set_id = context.object_sets->add_set(
        //            set_desc.name ); context.set_samplers.push_back(
        //                std::make_unique< geode::UniformPointSetSampler< 2 >
        //                >(
        //                    *context.domain ) );
        //            geode::add_birth_death_change_moves(
        //            context.set_samplers.back(),
        //                *proposal_kernel, set_id, set_desc.birth_ratio,
        //                set_desc.death_ratio, set_desc.change_ratio );
        //        }
        //        return proposal_kernel;

        for( const auto& set_cfg : config.sets )
        {
            const auto set_id = context.object_sets->add_set( set_cfg.name );
            geode_unused( set_id );
        }

        // -------------------------
        // Model
        // -------------------------
        context.model = geode::build_model< ObjectType >(
            config.model, *context.object_sets, *context.domain );

        // -------------------------
        //  Proposal
        // -------------------------
        auto proposal_kernel =
            std::make_unique< geode::ProposalKernel< ObjectType > >();
        for( const auto& set_proposal : config.proposals )
        {
            const auto set_id =
                context.object_sets->get_set_uuid( set_proposal.name );
            context.set_samplers.push_back(
                std::make_unique< geode::UniformPointSetSampler< 2 > >(
                    *context.domain ) );

            geode::add_birth_death_change_moves( context.set_samplers.back(),
                *proposal_kernel, set_id, set_proposal.birth_ratio,
                set_proposal.death_ratio, set_proposal.change_ratio );
        }

        // -------------------------
        // MH sampler
        // -------------------------
        context.mh_sampler =
            std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                *context.model, std::move( proposal_kernel ) );

        return context;
    }

} // namespace geode