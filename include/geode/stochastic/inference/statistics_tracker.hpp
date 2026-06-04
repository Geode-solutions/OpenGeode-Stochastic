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

#include <geode/basic/common.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/models/model.hpp>

namespace geode
{
    template < typename ObjectType >
    class StatisticsTracker
    {
    public:
        StatisticsTracker( const Model< ObjectType >& model ) : model_{ model }
        {
            means_.resize( model.nb_terms(), 0.0 );
            m2_.resize( model.nb_terms(), 0.0 );
        }

        [[nodiscard]] index_t statiscal_count() const
        {
            return count_;
        }

        void add_realization( const std::vector< double >& values )
        {
            ++count_;
            for( const auto value_id : geode::Range{ values.size() } )
            {
                auto& value = values[value_id];
                auto& mean = means_[value_id];
                auto& sum_of_squares = m2_[value_id];

                const auto delta = value - mean;
                mean += delta / count_;
                const auto delta2 = value - mean;
                sum_of_squares += delta * delta2;
            }
        }

        [[nodiscard]] double mean( const uuid& term_uuid ) const
        {
            return means_[model_.term_index( term_uuid )];
        }

        [[nodiscard]] const std::vector< double >& means() const
        {
            return means_;
        }

        [[nodiscard]] double variance( const uuid& term_uuid ) const
        {
            return variance( model_.term_index( term_uuid ) );
        }

        [[nodiscard]] std::vector< double > variances() const
        {
            std::vector< double > variances;
            variances.reserve( model_.nb_terms() );
            for( const auto variance_id : geode::Range{ model_.nb_terms() } )
            {
                variances.emplace_back( this->variance( variance_id ) );
            }
            return variances;
        }

    private:
        [[nodiscard]] double variance( index_t term_index ) const
        {
            if( count_ < 2 )
            {
                return 0.0;
            }
            return m2_[term_index] / ( count_ - 1 );
        }

    private:
        const Model< ObjectType >& model_;

        std::vector< double > means_;
        std::vector< double > m2_;
        index_t count_{ 0 };
    };
} // namespace geode