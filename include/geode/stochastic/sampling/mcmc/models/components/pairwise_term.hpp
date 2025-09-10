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
#include <geode/stochastic/spatial/object_set.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    // rename PairwiseInteractionTerm
    template < typename Type, typename InteractionFunc >
    class PairwiseTerm : public EnergyTerm< Type >
    {
    public:
        explicit PairwiseTerm( std::string_view name,
            double gamma,
            InteractionFunc interaction_func )
            : EnergyTerm< Type >( name, gamma ),
              interaction_func_( std::move( interaction_func ) )
        {
        }

        explicit PairwiseTerm( std::string_view name,
            double gamma,
            InteractionFunc interaction_func,
            const uuid& subset_id )
            : EnergyTerm< Type >( name, gamma, subset_id ),
              interaction_func_( std::move( interaction_func ) )
        {
        }

        double total_log( const ObjectSet< Type >& state ) const final
        {
            const auto interaction_weight = statistic( state );
            return this->contribution( interaction_weight );
        }

        double delta_log_add( const ObjectSet< Type >& state,
            const Type& new_object,
            const uuid& new_object_subset_id ) const final
        {
            if( !this->is_targeted_subset( new_object_subset_id ) )
            {
                return 0.0;
            }
            double delta = 0.0;
            const auto neighbors = state.neighbors( new_object, 1.1 );
            for( const auto& neigh_id : neighbors )
            {
                const auto& neigh_obj = state.get_object( neigh_id );
                delta += static_cast< double >( interaction_func_( neigh_obj,
                    neigh_id.subset, new_object, new_object_subset_id ) );
            }
            return this->contribution( delta );
        }

        double delta_log_remove(
            const ObjectSet< Type >& state, ObjectId object_id ) const override
        {
            if( !this->is_targeted_subset( object_id.subset ) )
            {
                return 0.0;
            }
            double delta = 0.0;
            const auto& to_remove = state.get_object( object_id );
            const auto neighbors = state.neighbors( object_id, 1.1 );
            for( auto neigh_id : neighbors )
            {
                const auto& neigh_obj = state.get_object( neigh_id );
                delta += static_cast< double >( interaction_func_(
                    neigh_obj, neigh_id.subset, to_remove, object_id.subset ) );
            }
            return this->contribution( -delta );
        }

        double delta_log_change( const ObjectSet< Type >& state,
            ObjectId old_object_id,
            const Type& new_object,
            const uuid& new_object_subset_id ) const override
        {
            if( !this->is_targeted_subset( old_object_id.subset )
                || !this->is_targeted_subset( new_object_subset_id ) )
            {
                return 0.0;
            }
            double delta = 0.0;

            // Remove old object's interactions
            const auto& old_obj = state.get_object( old_object_id );
            const auto old_neighbors = state.neighbors( old_object_id, 1.1 );
            for( auto neigh_id : old_neighbors )
            {
                const auto& neigh_obj = state.get_object( neigh_id );
                delta -= static_cast< double >( interaction_func_( neigh_obj,
                    neigh_id.subset, old_obj, old_object_id.subset ) );
            }

            // Add new object's interactions
            const auto new_neighbors = state.neighbors( new_object, 1.1 );
            for( auto neigh_id : new_neighbors )
            {
                if( old_object_id == neigh_id )
                    continue; // avoid double-counting
                const auto& neigh_obj = state.get_object( neigh_id );
                delta += static_cast< double >( interaction_func_( neigh_obj,
                    neigh_id.subset, new_object, new_object_subset_id ) );
            }

            return this->contribution( delta );
        }

        double statistic( const ObjectSet< Type >& state ) const override
        {
            double sum = 0.0;
            this->for_each_targeted_object( state, [&]( const ObjectId&
                                                           obj_id ) {
                const auto& cur_obj = state.get_object( obj_id );
                const auto neighbors =
                    state.get_all_object(); // state.neighbors( obj_id, 1.1 );
                for( const auto& neigh_obj_id : neighbors )
                {
                    if( neigh_obj_id == obj_id )
                    {
                        continue;
                    }
                    sum += static_cast< double >( interaction_func_( cur_obj,
                        obj_id.subset, state.get_object( neigh_obj_id ),
                        neigh_obj_id.subset ) );
                }
            } );
            return sum / 2.;
        }

    private:
        InteractionFunc interaction_func_;
    };
} // namespace geode