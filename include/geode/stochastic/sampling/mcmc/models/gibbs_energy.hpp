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

#include <geode/stochastic/sampling/mcmc/models/energy_term_collection.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace geode
{
    template < typename ObjectType >
    class GibbsEnergy
    {
    public:
        explicit GibbsEnergy(
            const EnergyTermCollection< ObjectType > &collection )
            : energy_terms_collection_( collection )
        {
        }

        double total_log_energy( const ObjectSets< ObjectType > &state ) const
        {
            double log_energy = 0.0;
            const auto &energy_terms = energy_terms_collection_.all_terms();
            for( auto &[id, term] : energy_terms )
            {
                log_energy += term->total_log( state );
            }
            return log_energy;
        }

        double total_log_energy_for_set(
            const ObjectSets< ObjectType > &state, const uuid &set_id ) const
        {
            double log_energy = 0.0;
            for( const auto term :
                energy_terms_collection_.terms_for_set( set_id ) )
            {
                log_energy += term->total_log( state );
            }
            return log_energy;
        }

        double delta_log_add( const ObjectSets< ObjectType > &state,
            const ObjectRef< ObjectType > &new_object ) const
        {
            double log_energy = 0.0;
            for( const auto &term :
                energy_terms_collection_.terms_for_set( new_object.set_id ) )
            {
                log_energy += term->delta_log_add( state, new_object );
            }
            return log_energy;
        }

        double delta_log_remove(
            const ObjectSets< ObjectType > &state, const ObjectId &id ) const
        {
            double log_energy = 0.0;
            for( const auto &term :
                energy_terms_collection_.terms_for_set( id.set_id ) )
            {
                log_energy += term->delta_log_remove( state, id );
            }
            return log_energy;
        }

        double delta_log_change( const ObjectSets< ObjectType > &state,
            const ObjectId &old_id,
            const ObjectRef< ObjectType > &new_object ) const
        {
            double log_energy = 0.0;
            for( const auto &term :
                energy_terms_collection_.terms_for_set( old_id.set_id ) )
            {
                log_energy +=
                    term->delta_log_change( state, old_id, new_object );
            }
            return log_energy;
        }

    private:
        const EnergyTermCollection< ObjectType > &energy_terms_collection_;
    };
} // namespace geode