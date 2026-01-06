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
#include <geode/stochastic/common.hpp>

namespace geode
{

    class StatisticsMonitor
    {
    public:
        StatisticsMonitor( StatisticsMonitor&& ) = default;
        StatisticsMonitor( const StatisticsMonitor& ) = default;
        StatisticsMonitor& operator=( StatisticsMonitor&& ) noexcept = default;
        StatisticsMonitor& operator=(
            const StatisticsMonitor& ) noexcept = default;

        StatisticsMonitor( const index_t nb_energy_terms )
        {
            means_.resize( nb_energy_terms, 0.0 );
            variances_.resize( nb_energy_terms, 0.0 );
        }

        void add_realization( const std::vector< double >& values )
        {
            OPENGEODE_EXCEPTION( values.size() == means_.size(),
                "[StatisticsMonitor] - Mismatch between realization size and "
                "expected number of statistics." );
            ++count_;
            for( size_t i = 0; i < values.size(); ++i )
            {
                double delta = values[i] - means_[i];
                means_[i] += delta / count_;
                if( count_ > 1 )
                    variances_[i] = ( ( count_ - 2 ) * variances_[i]
                                        + delta * ( values[i] - means_[i] ) )
                                    / ( count_ - 1 );
            }
        }

        const index_t statiscal_count() const
        {
            return count_;
        }

        const std::vector< double >& means() const
        {
            return means_;
        }

        const std::vector< double >& variances() const
        {
            return variances_;
        }

    private:
        std::vector< double > means_;
        std::vector< double > variances_;
        index_t count_{ 0 };
    };

} // namespace geode