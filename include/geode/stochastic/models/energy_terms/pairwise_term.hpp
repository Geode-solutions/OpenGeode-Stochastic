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

#include <absl/container/flat_hash_map.h>

#include <geode/basic/range.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>

#include <geode/stochastic/models/energy_terms/energy_term.hpp>

namespace geode
{
    template < typename ObjectType >
    class PairwiseTerm : public EnergyTerm< ObjectType >
    {
    public:
        explicit PairwiseTerm( std::string_view name,
            double gamma,
            std::vector< uuid >&& impacted_set_ids,
            const SpatialDomain< ObjectType::dim >& domain,
            absl::flat_hash_map< uuid, std::vector< uuid > >&&
                objectset_adjacency_map,
            std::unique_ptr< PairwiseInteraction< ObjectType > > interaction )
            : EnergyTerm< ObjectType >(
                  name, gamma, std::move( impacted_set_ids ), domain ),
              objectset_adjacency_map_( std::move( objectset_adjacency_map ) ),
              interaction_( std::move( interaction ) )
        {
        }

        [[nodiscard]] double total_log(
            const ObjectSets< ObjectType >& state ) const final
        {
            const auto interaction_weight = statistic( state );
            return this->contribution( interaction_weight );
        }

        [[nodiscard]] double delta_log_add(
            const ObjectSets< ObjectType >& state,
            const ObjectRef< ObjectType >& new_object ) const final
        {
            if( !this->is_impacted_set( new_object.set_id ) )
            {
                return 0.0;
            }
            auto delta = compute_local_interactions_with_neighbors(
                new_object, std::nullopt, state );
            return this->contribution( delta );
        }

        [[nodiscard]] double delta_log_remove(
            const ObjectSets< ObjectType >& state,
            const ObjectId& object_id ) const override
        {
            if( !this->is_impacted_set( object_id.set_id ) )
            {
                return 0.0;
            }
            ObjectRef< ObjectType > old_object{ state.get_object( object_id ),
                object_id.set_id };
            auto delta = compute_local_interactions_with_neighbors(
                old_object, object_id, state );
            return this->contribution( -delta );
        }

        [[nodiscard]] double delta_log_change(
            const ObjectSets< ObjectType >& state,
            const ObjectId& old_object_id,
            const ObjectRef< ObjectType >& new_object ) const override
        {
            auto delta_new = compute_local_interactions_with_neighbors(
                new_object, old_object_id, state );
            ObjectRef< ObjectType > old_object{
                state.get_object( old_object_id ), old_object_id.set_id
            };
            auto delta_old = compute_local_interactions_with_neighbors(
                old_object, old_object_id, state );
            auto delta = delta_new - delta_old;
            return this->contribution( delta );
        }

        [[nodiscard]] double statistic(
            const ObjectSets< ObjectType >& state ) const override
        {
            double sum = 0.0;
            this->for_each_object_in_sets( state, this->impacted_set_ids(),
                [&sum, &state, this]( const ObjectId& cur_obj_id ) {
                    sum += this->accumulate_interactions_with_neighbors(
                        cur_obj_id, state );
                } );
            return sum;
        }

    private:
        double compute_local_interactions_with_neighbors(
            const ObjectRef< ObjectType >& object_ref,
            std::optional< ObjectId > exclude_id,
            const ObjectSets< ObjectType >& state ) const
        {
            const auto impacted_set_it =
                objectset_adjacency_map_.find( object_ref.set_id );
            if( impacted_set_it == objectset_adjacency_map_.end() )
            {
                return 0.;
            }
            const auto neighbors = state.neighbors( object_ref.object,
                impacted_set_it->second,
                interaction_->neighborhood_searching_distance(), exclude_id );
            double sum = 0.0;
            for( const auto& neigh_id : neighbors )
            {
                ObjectRef< ObjectType > neigh_object{
                    state.get_object( neigh_id ), neigh_id.set_id
                };
                if( !is_any_in_domain(
                        object_ref.object, neigh_object.object ) )
                {
                    continue;
                }
                sum += interaction_->evaluate( object_ref, neigh_object );
            }
            return sum;
        }

        double accumulate_interactions_with_neighbors(
            const ObjectId& object_id,
            const ObjectSets< ObjectType >& state ) const
        {
            const auto& cur_obj = state.get_object( object_id );
            if( !is_in_domain( cur_obj ) )
            {
                return 0.;
            }
            const auto impacted_set_it =
                objectset_adjacency_map_.find( object_id.set_id );
            if( impacted_set_it == objectset_adjacency_map_.end() )
            {
                return 0.;
            }
            const auto neighbors = state.neighbors( cur_obj,
                impacted_set_it->second,
                interaction_->neighborhood_searching_distance(), object_id );
            ObjectRef< ObjectType > object_ref{ cur_obj, object_id.set_id };
            double sum = 0.0;
            for( const auto& neigh_id : neighbors )
            {
                const auto& neigh_obj = state.get_object( neigh_id );
                if( neigh_id < object_id && is_in_domain( neigh_obj ) )
                {
                    continue;
                }
                ObjectRef< ObjectType > neigh_object{ neigh_obj,
                    neigh_id.set_id };
                sum += interaction_->evaluate( object_ref, neigh_object );
            }
            return sum;
        }

        bool is_any_in_domain(
            const ObjectType& object1, const ObjectType& object2 ) const
        {
            return is_in_domain( object1 ) || is_in_domain( object2 );
        }

        bool is_in_domain( const ObjectType& object ) const
        {
            return SpatialDomainChecker< ObjectType >::is_anchored_in_domain(
                this->domain(), object );
        }

    private:
        absl::flat_hash_map< uuid, std::vector< uuid > >
            objectset_adjacency_map_;
        std::unique_ptr< PairwiseInteraction< ObjectType > > interaction_;
    };
} // namespace geode