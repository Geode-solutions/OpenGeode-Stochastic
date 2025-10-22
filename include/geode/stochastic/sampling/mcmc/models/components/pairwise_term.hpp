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

#include <geode/basic/range.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    // rename PairwiseInteractionTerm
    template < typename ObjectType >
    class PairwiseTerm : public EnergyTerm< ObjectType >
    {
    public:
        explicit PairwiseTerm( std::string_view name,
            double gamma,
            absl::flat_hash_set< uuid > targeted_set_ids,
            std::unique_ptr< PairwiseInteraction< ObjectType > > interaction )
            : EnergyTerm< ObjectType >(
                  name, gamma, std::move( targeted_set_ids ) ),
              interaction_( std::move( interaction ) )
        {
        }

        double total_log( const ObjectSets< ObjectType >& state ) const final
        {
            const auto interaction_weight = statistic( state );
            return this->contribution( interaction_weight );
        }

        double delta_log_add( const ObjectSets< ObjectType >& state,
            const ObjectRef< ObjectType >& new_object ) const final
        {
            if( !this->is_targeted_set( new_object.set_id ) )
            {
                return 0.0;
            }
            double delta = 0.0;
            const auto neighbors = state.neighbors( new_object.object,
                interaction_->neighborhood_searching_distance() );
            for( const auto& neigh_id : neighbors )
            {
                geode::ObjectRef< ObjectType > neigh_object{
                    state.get_object( neigh_id ), neigh_id.set_id
                };

                delta += interaction_->evaluate( new_object, neigh_object );
            }
            return this->contribution( delta );
        }

        double delta_log_remove( const ObjectSets< ObjectType >& state,
            const ObjectId& object_id ) const override
        {
            if( !this->is_targeted_set( object_id.set_id ) )
            {
                return 0.0;
            }
            double delta = 0.0;
            ObjectRef< ObjectType > object_to_remove{
                state.get_object( object_id ), object_id.set_id
            };
            const auto neighbors = state.neighbors(
                object_id, interaction_->neighborhood_searching_distance() );
            for( auto neigh_id : neighbors )
            {
                ObjectRef< ObjectType > neigh_object{
                    state.get_object( neigh_id ), neigh_id.set_id
                };
                delta +=
                    interaction_->evaluate( object_to_remove, neigh_object );
            }
            return this->contribution( -delta );
        }

        double delta_log_change( const ObjectSets< ObjectType >& state,
            const ObjectId& old_object_id,
            const ObjectRef< ObjectType >& new_object ) const override
        {
            if( !this->is_targeted_set( old_object_id.set_id )
                || !this->is_targeted_set( new_object.set_id ) )
            {
                return 0.0;
            }
            double delta = 0.0;

            // Remove old object's interactions
            ObjectRef< ObjectType > object_to_remove{
                state.get_object( old_object_id ), old_object_id.set_id
            };
            const auto old_neighbors = state.neighbors( old_object_id,
                interaction_->neighborhood_searching_distance() );
            for( auto neigh_id : old_neighbors )
            {
                ObjectRef< ObjectType > neigh_object{
                    state.get_object( neigh_id ), neigh_id.set_id
                };
                delta -=
                    interaction_->evaluate( object_to_remove, neigh_object );
            }

            // Add new object's interactions
            const auto new_neighbors = state.neighbors( new_object.object,
                interaction_->neighborhood_searching_distance() );
            for( auto neigh_id : new_neighbors )
            {
                if( old_object_id == neigh_id )
                {
                    continue; // avoid double-counting
                }
                ObjectRef< ObjectType > neigh_object{
                    state.get_object( neigh_id ), neigh_id.set_id
                };
                delta += interaction_->evaluate( new_object, neigh_object );
            }

            return this->contribution( delta );
        }

        double statistic( const ObjectSets< ObjectType >& state ) const override
        {
            double sum = 0.0;
            this->for_each_targeted_object( state, [&]( const ObjectId&
                                                           obj_id ) {
                const auto& cur_obj = state.get_object( obj_id );
                const auto neighbors =
                    state.get_all_object(); // state.neighbors( obj_id, 1.1 );
                for( const auto& neigh_obj_id : neighbors )
                {
                    // if( neigh_obj_id.set_id < obj_id.set_id )
                    //{
                    //     continue;
                    // }
                    // if( neigh_obj_id.set_id == obj_id.set_id
                    //     && neigh_obj_id.index <= obj_id.index )
                    //{
                    //     continue;
                    // }

                    if( neigh_obj_id == obj_id )
                    {
                        continue;
                    }
                    ObjectRef< ObjectType > object{ cur_obj, obj_id.set_id };
                    ObjectRef< ObjectType > neigh_object{
                        state.get_object( neigh_obj_id ), neigh_obj_id.set_id
                    };

                    sum += interaction_->evaluate( object, neigh_object );
                }
            } );
            return sum / 2.;
        }

    private:
        std::unique_ptr< PairwiseInteraction< ObjectType > > interaction_;
    };
} // namespace geode