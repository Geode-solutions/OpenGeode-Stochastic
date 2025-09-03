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

#include <geode/stochastic/configuration/configuration.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    template < typename Object >
    class IntensityTerm : public EnergyTerm< Object >
    {
    public:
        explicit IntensityTerm( double lambda, GroupId group_id )
            : EnergyTerm< Object >( lambda ), group_id_{ group_id }
        {
        }

        double total_log( const Configuration< Object >& state ) const final
        {
            const auto n = static_cast< double >( number_of_objects( state ) );
            return this->neg_log_parameter_.scale( n );
        }

        double delta_log_add( const Configuration< Object >& state,
            const Object& sample,
            const GroupId target_group_id ) const final
        {
            return target_group_id == group_id_
                       ? this->neg_log_parameter_.scale( 1. )
                       : 0.;
        }

        double delta_log_remove( const Configuration< Object >& state,
            ObjectId object_id ) const final
        {
            return object_id.group == group_id_
                       ? this->neg_log_parameter_.scale( -1. )
                       : 0.;
        }

        double delta_log_change( const Configuration< Object >& state,
            ObjectId old_object_id,
            const Object& new_sample ) const final
        {
            return 0.0;
        }

        double statistic( const Configuration< Object >& state ) const final
        {
            return static_cast< double >( number_of_objects( state ) );
        }

    private:
        index_t number_of_objects( const Configuration< Object >& state ) const
        {
            return state.nb_objects_in_group( group_id_ );
        }

    private:
        GroupId group_id_;
    };
} // namespace geode