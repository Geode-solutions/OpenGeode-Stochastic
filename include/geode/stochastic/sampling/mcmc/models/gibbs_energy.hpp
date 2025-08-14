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
    template < typename Geometry >
    class GibbsEnergy
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( GibbsEnergy );

    public:
        GibbsEnergy() = default;
        ~GibbsEnergy() = default;

        void add_energy_term( std::unique_ptr< EnergyTerm< Geometry > > term )
        {
            energy_terms_.push_back( std::move( term ) );
        }
        void clear_energy_terms()
        {
            energy_terms_.clear();
            energy_terms_.shrink_to_fit();
        }
        index_t number_of_energy_terms() const
        {
            return energy_terms_.size();
        }
        std::vector< double > ordered_energy_term_parameter() const
        {
            std::vector< double > values;
            values.reserve( energy_terms_.size() );
            for( const auto& term : energy_terms_ )
            {
                values.emplace_back( term->parameter() );
            }
            return values;
        }

        std::string ordered_energy_term_parameter_string() const
        {
            const auto values = ordered_energy_term_parameter();
            return vector_string( values );
        }

        double log_energy_total( const Configuration< Geometry > state )
        {
            double log_energy{ 0.0 };
            for( const auto& term : energy_terms_ )
            {
                log_energy += term->log_total( state );
            }
            return log_energy;
        }

        double log_energy_delta_add( const Configuration< Geometry > state,
            const MarkedObject< Geometry >& sample )
        {
            double log_energy{ 0.0 };
            for( const auto& term : energy_terms_ )
            {
                log_energy += term->log_delta_add( state, sample );
            }
            return log_energy;
        }

        double log_energy_delta_remove(
            const Configuration< Geometry > state, index_t sample_id )
        {
            double log_energy{ 0.0 };
            for( const auto& term : energy_terms_ )
            {
                log_energy += term->log_delta_remove( state, sample_id );
            }
            return log_energy;
        }

        double log_energy_delta_change( const Configuration< Geometry > state,
            index_t old_sample_id,
            const MarkedObject< Geometry >& new_sample )
        {
            double log_energy{ 0.0 };
            for( const auto& term : energy_terms_ )
            {
                log_energy +=
                    term->log_delta_change( state, old_sample_id, new_sample );
            }
            return log_energy;
        }

        std::vector< double > ordered_energy_term_statistics(
            const Configuration< Geometry > state ) const
        {
            std::vector< double > values;
            values.reserve( energy_terms_.size() );
            for( const auto& term : energy_terms_ )
            {
                values.emplace_back(
                    static_cast< double >( term->statistic( state ) ) );
            }
            return values;
        }

        std::string ordered_energy_term_statistics_string(
            const Configuration< Geometry > state ) const
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
        std::vector< std::unique_ptr< EnergyTerm< Geometry > > > energy_terms_;
    };

} // namespace geode