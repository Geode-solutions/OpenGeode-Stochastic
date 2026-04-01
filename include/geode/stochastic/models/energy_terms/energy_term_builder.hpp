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
#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_density_term(
        const DensityTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        auto set_ids =
            object_sets.get_existing_set_uuids( cfg.object_set_names );

        return std::make_unique< geode::DensityTerm< ObjectType > >(
            cfg.term_name, cfg.lambda, set_ids, domain );
    }

    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_pairwise_term(
        const PairwiseTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        auto set_ids =
            object_sets.get_existing_set_uuids( cfg.object_set_names );

        auto interaction =
            std::make_unique< geode::MinimalDistanceCutoff< ObjectType > >( 0.5
                /*,interaction_desc.interaction_scope*/ );

        return std::make_unique< geode::PairwiseTerm< geode::Point2D > >(
            cfg.term_name, cfg.gamma, set_ids, std::move( interaction ),
            domain );
    }

    template < typename ObjectType >
    std::unique_ptr< EnergyTerm< ObjectType > > build_term(
        const EnergyTermConfig& cfg,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        return std::visit(
            [&]( auto&& c ) -> std::unique_ptr< EnergyTerm< ObjectType > > {
                using T = std::decay_t< decltype( c ) >;

                if constexpr( std::is_same_v< T, DensityTermConfig > )
                    return build_density_term< ObjectType >(
                        c, object_sets, domain );

                else if constexpr( std::is_same_v< T, PairwiseTermConfig > )
                    return build_pairwise_term< ObjectType >(
                        c, object_sets, domain );
            },
            cfg );
    }

} // namespace geode