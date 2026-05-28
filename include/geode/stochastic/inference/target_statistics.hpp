/*
 * Copyright (c) 2019 - 2026 Geode-solutions
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
#include <geode/basic/uuid.hpp>

namespace geode
{
    struct TargetStatisticConfig
    {
        std::string term_name;
        double value;
        double tolerance;
    };

    template < typename ObjectType >
    class TargetStatistics
    {
    public:
        explicit TargetStatistics( const Model< ObjectType >& model,
            const std::vector< TargetStatisticConfig >& statistic_targets )
            : model_( model )
        {
            values_.resize( model.nb_terms(), 0.0 );
            tolerances_.resize( model.nb_terms(), 0.0 );
            active_.resize( model.nb_terms(), false );
            for( const auto& target : statistic_targets )
            {
                set_target( target );
            }
        }

        const Model< ObjectType >& model() const
        {
            return model_;
        }

        bool has_target( const uuid& term_uuid ) const
        {
            return active_[model_.term_index( term_uuid )];
        }

        double target( const uuid& term_uuid ) const
        {
            return values_[model_.term_index( term_uuid )];
        }

        double tolerance( const uuid& term_uuid ) const
        {
            return tolerances_[model_.term_index( term_uuid )];
        }

        std::vector< uuid > active_terms() const
        {
            std::vector< uuid > active_terms_uuid;

            for( const auto& term : model_.terms().energy_terms() )
            {
                const auto& id = term->id();
                if( active_[model_.term_index( id )] )
                {
                    active_terms_uuid.push_back( id );
                }
            }
            return active_terms_uuid;
        }

    private:
        void set_target( const TargetStatisticConfig& statistic )
        {
            const auto term_uuid =
                model_.terms().get_term_uuid( statistic.term_name );
            const auto idx = model_.term_index( term_uuid );

            values_[idx] = statistic.value;
            tolerances_[idx] = statistic.tolerance;
            active_[idx] = true;
        }

    private:
        const Model< ObjectType >& model_;

        std::vector< double > values_;
        std::vector< double > tolerances_;
        std::vector< bool > active_;
    };
} // namespace geode
