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
    template < typename Geometry, typename InteractionFunc >
    class PairwiseTerm : public EnergyTerm< Geometry >
    {
    public:
        explicit PairwiseTerm( double gamma, InteractionFunc func )
            : EnergyTerm< Geometry >( gamma ),
              interaction_func_( std::move( func ) )
        {
        }

        double total_log( const Configuration< Geometry >& state ) const final
        {
            double interaction_weight = 0.0;
            for( const auto obj_id : geode::Range{ state.size() } )
            {
                for( const auto neigh_id :
                    geode::Range{ obj_id + 1, state.size() } )
                {
                    interaction_weight += static_cast< double >(
                        interaction_func_( state[obj_id], state[neigh_id] ) );
                }
            }
            return this->neg_log_parameter_.scale( interaction_weight );
        }

        double delta_log_add( const Configuration< Geometry >& state,
            const MarkedObject< Geometry >& new_object ) const final
        {
            double interaction_weight = 0.0;
            for( const auto obj_id : geode::Range( state.size() ) )
            {
                interaction_weight += static_cast< double >(
                    interaction_func_( state[obj_id], new_object ) );
            }
            return this->neg_log_parameter_.scale( interaction_weight );
        }

        double delta_log_remove( const Configuration< Geometry >& state,
            index_t object_id ) const final
        {
            const auto& removed = state[object_id];
            double interaction_weight = 0.0;
            for( const auto neigh_id : geode::Range{ state.size() } )
            {
                if( neigh_id == object_id )
                {
                    continue;
                }
                interaction_weight += static_cast< double >(
                    interaction_func_( removed, state[neigh_id] ) );
            }
            return this->neg_log_parameter_.scale( -interaction_weight );
        }

        double delta_log_change( const Configuration< Geometry >& state,
            index_t old_object_id,
            const MarkedObject< Geometry >& new_object ) const final
        {
            const auto& old_obj = state[old_object_id];
            double interaction_weight_remove = 0.0;
            double interaction_weight_add = 0.0;

            for( const auto neigh_id : geode::Range{ state.size() } )
            {
                if( neigh_id == old_object_id )
                {
                    continue;
                }
                interaction_weight_remove += static_cast< double >(
                    interaction_func_( old_obj, state[neigh_id] ) );
                interaction_weight_add += static_cast< double >(
                    interaction_func_( new_object, state[neigh_id] ) );
            }
            return this->neg_log_parameter_.scale(
                interaction_weight_add - interaction_weight_remove );
        }

        double statistic( const Configuration< Geometry >& state ) const final
        {
            double sum = 0.0;
            for( const auto obj_id : geode::Range{ state.size() } )
            {
                for( const auto neigh_id :
                    geode::Range{ obj_id + 1, state.size() } )
                {
                    sum += interaction_func_( state[obj_id], state[neigh_id] );
                }
            }
            return sum;
        }

    private:
        InteractionFunc interaction_func_;
    };
} // namespace geode