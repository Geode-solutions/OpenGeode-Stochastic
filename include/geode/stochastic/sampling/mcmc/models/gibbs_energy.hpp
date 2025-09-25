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

#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>
#include <geode/stochastic/spatial/object_set.hpp>

#include <absl/container/flat_hash_map.h>
namespace geode
{
    template < typename Type >
    class GibbsEnergy
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( GibbsEnergy );

    public:
        GibbsEnergy() = default;
        ~GibbsEnergy() = default;

        void add_energy_term( std::unique_ptr< EnergyTerm< Type > > term )
        {
            energy_terms_.emplace( term->id(), std::move( term ) );
        }
        void clear_energy_terms()
        {
            energy_terms_.clear();
        }
        index_t number_of_energy_terms() const
        {
            return energy_terms_.size();
        }
        std::vector< double > ordered_energy_term_parameter() const
        {
            std::vector< double > values;
            values.reserve( energy_terms_.size() );
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                values.emplace_back( term->parameter() );
            }
            return values;
        }

        std::string ordered_energy_term_parameter_string() const
        {
            const auto values = ordered_energy_term_parameter();
            return vector_string( values );
        }

        double total_log_energy( const ObjectSet< Type >& state ) const
        {
            double log_energy{ 0.0 };
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                log_energy += term->total_log( state );
            }
            return log_energy;
        }

        double delta_log_energy_add( const ObjectSet< Type >& state,
            const ObjectRef< Type >& new_object ) const
        {
            double log_energy{ 0.0 };
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                log_energy += term->delta_log_add( state, new_object );
            }
            return log_energy;
        }

        double delta_log_energy_remove(
            const ObjectSet< Type >& state, ObjectId object_id ) const
        {
            double log_energy{ 0.0 };
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                log_energy += term->delta_log_remove( state, object_id );
            }
            return log_energy;
        }

        double delta_log_energy_change( const ObjectSet< Type >& state,
            ObjectId old_object_id,
            const Type& new_object,
            const uuid& new_object_subset_id ) const
        {
            double log_energy{ 0.0 };
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                log_energy += term->delta_log_change(
                    state, old_object_id, new_object, new_object_subset_id );
            }
            return log_energy;
        }

        std::vector< double > ordered_energy_term_statistics(
            const ObjectSet< Type >& state ) const
        {
            std::vector< double > values;
            values.reserve( energy_terms_.size() );
            for( const auto& [id, term] : energy_terms_ )
            {
                geode_unused( id );
                values.emplace_back(
                    static_cast< double >( term->statistic( state ) ) );
            }
            return values;
        }

        std::string ordered_energy_term_statistics_string(
            const ObjectSet< Type >& state ) const
        {
            const auto stats = ordered_energy_term_statistics( state );
            return vector_string( stats );
        }

    private:
        std::string vector_string( const std::vector< double >& vector ) const
        {
            std::string str_values;
            for( const auto& value : vector )
            {
                absl::StrAppend( &str_values, " ", value );
            }
            absl::StrAppend( &str_values, "\n" );
            return str_values;
        }

    private:
        absl::flat_hash_map< geode::uuid,
            std::unique_ptr< EnergyTerm< Type > > >
            energy_terms_;
    };

} // namespace geode