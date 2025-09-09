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

#include <geode/stochastic/spatial/object_set.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    template < typename Type >
    class IntensityTerm : public EnergyTerm< Type >
    {
    public:
        explicit IntensityTerm( double lambda, uuid subset_id )
            : EnergyTerm< Type >( lambda ), subset_id_{ subset_id }
        {
        }

        double total_log( const ObjectSet< Type >& state ) const final
        {
            const auto n = static_cast< double >( number_of_objects( state ) );
            return this->neg_log_parameter_.scale( n );
        }

        double delta_log_add( const ObjectSet< Type >& state,
            const Type& sample,
            const uuid target_subset_id ) const final
        {
            return target_subset_id == subset_id_
                       ? this->neg_log_parameter_.scale( 1. )
                       : 0.;
        }

        double delta_log_remove(
            const ObjectSet< Type >& state, ObjectId object_id ) const final
        {
            return object_id.subset == subset_id_
                       ? this->neg_log_parameter_.scale( -1. )
                       : 0.;
        }

        double delta_log_change( const ObjectSet< Type >& state,
            ObjectId old_object_id,
            const Type& new_sample ) const final
        {
            return 0.0;
        }

        double statistic( const ObjectSet< Type >& state ) const final
        {
            return static_cast< double >( number_of_objects( state ) );
        }

    private:
        index_t number_of_objects( const ObjectSet< Type >& state ) const
        {
            return state.nb_objects_in_subset( subset_id_ );
        }

    private:
        uuid subset_id_;
    };
} // namespace geode