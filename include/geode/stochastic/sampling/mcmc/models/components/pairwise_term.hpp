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
#include <geode/stochastic/configuration/configuration.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    template < typename Type, typename InteractionFunc >
    class PairwiseTerm : public EnergyTerm< Type >
    {
    public:
        explicit PairwiseTerm( double gamma, InteractionFunc func )
            : EnergyTerm< Type >( gamma ),
              interaction_func_( std::move( func ) )
        {
        }

        double total_log( const ObjectSet< Type >& state ) const final
        {
            const auto interaction_weight = statistic( state );
            return this->neg_log_parameter_.scale( interaction_weight );
        }

        double delta_log_add( const ObjectSet< Type >& state,
            const Type& new_object,
            uuid subset_id ) const final
        {
            geode_unused( subset_id );
            const auto neighbors = state.neighbors( new_object, 1.1 );
            double interaction_weight = 0.0;
            for( const auto neigh_obj_id : neighbors )
            {
                interaction_weight += static_cast< double >( interaction_func_(
                    state.get_object( neigh_obj_id ), new_object ) );
            }
            return this->neg_log_parameter_.scale( interaction_weight );
        }

        double delta_log_remove(
            const ObjectSet< Type >& state, ObjectId object_id ) const final
        {
            const auto& to_remove = state.get_object( object_id );
            const auto neighbors = state.neighbors( object_id, 1.1 );
            double interaction_weight = 0.0;
            for( const auto neigh_obj_id : neighbors )
            {
                interaction_weight += static_cast< double >( interaction_func_(
                    to_remove, state.get_object( neigh_obj_id ) ) );
            }
            return this->neg_log_parameter_.scale( -interaction_weight );
        }

        double delta_log_change( const ObjectSet< Type >& state,
            ObjectId old_object_id,
            const Type& new_object ) const final
        {
            const auto new_neighbors = state.neighbors( new_object, 1.1 );
            double interaction_weight_add = 0.0;
            for( const auto neigh_obj_id : new_neighbors )
            {
                if( old_object_id == neigh_obj_id )
                {
                    continue;
                }
                interaction_weight_add +=
                    static_cast< double >( interaction_func_(
                        state.get_object( neigh_obj_id ), new_object ) );
            }

            const auto& to_remove = state.get_object( old_object_id );
            const auto old_neighbors = state.neighbors( old_object_id, 1.1 );
            double interaction_weight_remove = 0.0;
            for( const auto neigh_obj_id : old_neighbors )
            {
                interaction_weight_remove +=
                    static_cast< double >( interaction_func_(
                        to_remove, state.get_object( neigh_obj_id ) ) );
            }

            return this->neg_log_parameter_.scale(
                interaction_weight_add - interaction_weight_remove );
        }

        double statistic( const ObjectSet< Type >& state ) const final
        {
            double sum = 0.0;
            const auto all_object_ids = state.get_all_object();
            for( const auto obj_id : geode::Range{ all_object_ids.size() } )
            {
                const auto& cur_objet =
                    state.get_object( all_object_ids[obj_id] );
                const auto neighbors =
                    state.neighbors( all_object_ids[obj_id], 1.1 );
                for( const auto neigh_obj_id :
                    geode::Range{ neighbors.size() } )
                {
                    if( all_object_ids[obj_id].index
                        > neighbors[neigh_obj_id].index )
                    {
                        sum +=
                            static_cast< double >( interaction_func_( cur_objet,
                                state.get_object( neighbors[neigh_obj_id] ) ) );
                    }
                }
            }
            return sum;
        }

    private:
        InteractionFunc interaction_func_;
    };
} // namespace geode