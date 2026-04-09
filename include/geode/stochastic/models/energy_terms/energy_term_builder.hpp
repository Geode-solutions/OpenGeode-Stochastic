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

#include <absl/strings/str_join.h>

#include <geode/basic/uuid.hpp>

#include <string>
#include <variant>
#include <vector>

#include <geode/stochastic/models/energy_terms/density_term.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>
#include <geode/stochastic/models/energy_terms/pairwise_term.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_builder.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>

#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions_builder.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions_config.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_single_term(
        const SingleObjectTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        auto set_ids = object_sets.get_set_uuids( cfg.object_set_names );

        auto object_feature =
            build_single_object_feature< ObjectType >( cfg.object_feature );

        return std::make_unique< geode::SingleObjectTerm< ObjectType > >(
            cfg.term_name, cfg.lambda, std::move( set_ids ), domain,
            std::move( object_feature ) );
    }

    std::pair< std::vector< geode::uuid >,
        absl::flat_hash_map< uuid, std::vector< uuid > > >
        pairwise_builder_initialize_interactions_helper(
            const std::vector< std::pair< uuid, uuid > >& interacting_sets )
    {
        std::vector< geode::uuid > interacting_set_ids;
        absl::flat_hash_map< uuid, std::vector< uuid > >
            objectset_adjacency_map;
        interacting_set_ids.reserve( 2 * interacting_sets.size() );
        objectset_adjacency_map.reserve( 2 * interacting_sets.size() );
        for( const auto& [set1, set2] : interacting_sets )
        {
            objectset_adjacency_map[set1].push_back( set2 );
            objectset_adjacency_map[set2].push_back( set1 );

            interacting_set_ids.push_back( set1 );
            interacting_set_ids.push_back( set2 );
        }

        std::sort( interacting_set_ids.begin(), interacting_set_ids.end() );
        interacting_set_ids.erase( std::unique( interacting_set_ids.begin(),
                                       interacting_set_ids.end() ),
            interacting_set_ids.end() );
        interacting_set_ids.shrink_to_fit();

        for( auto& [_, adjacent_set_uuids] : objectset_adjacency_map )
        {
            geode_unused( _ );
            std::sort( adjacent_set_uuids.begin(), adjacent_set_uuids.end() );
            adjacent_set_uuids.erase( std::unique( adjacent_set_uuids.begin(),
                                          adjacent_set_uuids.end() ),
                adjacent_set_uuids.end() );
            adjacent_set_uuids.shrink_to_fit();
        }
        return std::make_pair( interacting_set_ids, objectset_adjacency_map );
    }
    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_pairwise_term(
        const PairwiseTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        auto set_id_interactions =
            object_sets.get_set_uuid_pairs( cfg.object_set_names_interactions );

        auto [interacting_set_ids, adjacent_set_uuids] =
            pairwise_builder_initialize_interactions_helper(
                set_id_interactions );

        auto interaction =
            build_pairwise_interaction< ObjectType >( cfg.interaction_config );

        return std::make_unique< geode::PairwiseTerm< ObjectType > >(
            cfg.term_name, cfg.gamma, std::move( interacting_set_ids ), domain,
            std::move( adjacent_set_uuids ), std::move( interaction ) );
    }

    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_energy_term_impl(
        const SingleObjectTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        return build_single_term< ObjectType >( cfg, object_sets, domain );
    }

    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_energy_term_impl(
        const PairwiseTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        return build_pairwise_term< ObjectType >( cfg, object_sets, domain );
    }

    template < typename ObjectType, typename NewEnergyTypeConfig >
    std::unique_ptr< EnergyTerm< ObjectType > > build_energy_term_impl(
        const NewEnergyTypeConfig&,
        const ObjectSets< ObjectType >&,
        const SpatialDomain< ObjectType::dim >& )
    {
        static_assert( sizeof( NewEnergyTypeConfig ) == 0,
            "Unsupported EnergyTermConfig type" );
        return nullptr;
    }
    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_energy_term_impl(
        const std::monostate&,
        const ObjectSets< ObjectType >&,
        const SpatialDomain< ObjectType::dim >& )
    {
        throw OpenGeodeException(
            "[EnergyTermBuilder] energy term config not initialized" );
    }
    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_energy_term(
        const EnergyTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        return std::visit(
            [&]( const auto& term_cfg )
                -> std::unique_ptr< EnergyTerm< ObjectType > > {
                return build_energy_term_impl< ObjectType >(
                    term_cfg, object_sets, domain );
            },
            cfg );
    }

} // namespace geode