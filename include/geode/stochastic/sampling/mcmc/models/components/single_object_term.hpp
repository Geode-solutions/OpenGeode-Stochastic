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

#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < typename ObjectType, typename ObjectContributionFunc >
    class SingleObjectTerm : public EnergyTerm< ObjectType >
    {
    public:
        explicit SingleObjectTerm( std::string_view name,
            double lambda,
            std::vector< uuid > targeted_set_ids,
            double scale,
            ObjectContributionFunc contribution_func,
            const SpatialDomain< ObjectType::dim >& domain )
            : EnergyTerm< ObjectType >(
                  name, lambda, std::move( targeted_set_ids ), domain ),
              scale_( scale ),
              contribution_func_( std::move( contribution_func ) )
        {
        }

        double total_log( const ObjectSets< ObjectType >& state ) const override
        {
            return this->contribution( scale_ * statistic( state ) );
        }

        double delta_log_add( const ObjectSets< ObjectType >& /*state*/,
            const ObjectRef< ObjectType >& new_object ) const override
        {
            if( !this->is_targeted_set( new_object.set_id ) )
            {
                return 0.0;
            }
            return this->contribution(
                scale_
                * contribution_func_( new_object.object, this->domain() ) );
        }

        double delta_log_remove( const ObjectSets< ObjectType >& state,
            const ObjectId& object_id ) const override
        {
            if( !this->is_targeted_set( object_id.set_id ) )
            {
                return 0.0;
            }
            return this->contribution(
                -scale_
                * contribution_func_(
                    state.get_object( object_id ), this->domain() ) );
        }

        double delta_log_change( const ObjectSets< ObjectType >& state,
            const ObjectId& old_object_id,
            const ObjectRef< ObjectType >& new_object ) const override
        {
            double delta =
                contribution_func_( new_object.object, this->domain() )
                - contribution_func_(
                    state.get_object( old_object_id ), this->domain() );
            return this->contribution( scale_ * delta );
        }

        double statistic( const ObjectSets< ObjectType >& state ) const override
        {
            double sum = 0.0;
            this->for_each_targeted_object(
                state, [&]( const ObjectId& obj_id ) {
                    const auto& obj = state.get_object( obj_id );
                    sum += contribution_func_( obj, this->domain() );
                } );
            return sum;
        }

    private:
        double scale_{ 1. };
        ObjectContributionFunc contribution_func_;
    };
} // namespace geode