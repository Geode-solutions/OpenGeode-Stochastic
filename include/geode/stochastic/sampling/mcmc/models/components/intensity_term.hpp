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
    /// ObjectCountTerm
    template < typename Type >
    class IntensityTerm : public EnergyTerm< Type >
    {
    public:
        explicit IntensityTerm( std::string_view name, double lambda )
            : EnergyTerm< Type >( name, lambda )
        {
        }

        explicit IntensityTerm(
            std::string_view name, double lambda, const uuid& subset_id )
            : EnergyTerm< Type >( name, lambda, subset_id )
        {
        }

        double total_log( const ObjectSet< Type >& state ) const override
        {
            const auto n = this->statistic( state );
            return this->contribution( n );
        }

        double delta_log_add( const ObjectSet< Type >& /*state*/,
            const ObjectRef< Type >& new_object ) const override
        {
            if( !this->is_targeted_subset( new_object.subset ) )
            {
                return 0.0;
            }
            return this->contribution( 1.0 );
        }

        double delta_log_remove( const ObjectSet< Type >& /*state*/,
            const ObjectId& object_id ) const override
        {
            if( !this->is_targeted_subset( object_id.subset ) )
            {
                return 0.0;
            }
            return this->contribution( -1.0 );
        }

        double delta_log_change( const ObjectSet< Type >& /*state*/,
            const ObjectId& /*old_object_id*/,
            const ObjectRef< Type >& /*new_object*/ ) const override
        {
            return 0.0;
        }

        double statistic( const ObjectSet< Type >& state ) const override
        {
            if( this->targeted_subset_id() )
            {
                return static_cast< double >( state.nb_objects_in_subset(
                    this->targeted_subset_id().value() ) );
            }
            return static_cast< double >( state.nb_objects() );
        }
    };
} // namespace geode