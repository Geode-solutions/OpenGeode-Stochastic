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
#include <geode/stochastic/common.hpp>
// #include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
// #include <geode/stochastic/sampling/mcmc/models/energy_term_collection.hpp>

#include <absl/strings/str_join.h>
#include <filesystem>
#include <fstream>

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
            sum.resize( nb_energy_terms, 0.0 );
            sum_squares.resize( nb_energy_terms, 0.0 );
            means.resize( nb_energy_terms, 0.0 );
            variances.resize( nb_energy_terms, 0.0 );
        }

        void add_realization( const std::vector< double >& values )
        {
            for( const auto stat_id : Range{ values.size() } )
            {
                sum[stat_id] += values[stat_id];
                sum_squares[stat_id] += values[stat_id] * values[stat_id];
            }
        }

        void finalize( const index_t nb_realizations )
        {
            for( const auto stat_id : Range{ sum.size() } )
            {
                means[stat_id] = sum[stat_id] / nb_realizations;
                double variance =
                    ( sum_squares[stat_id]
                        - ( sum[stat_id] * sum[stat_id] ) / nb_realizations )
                    / ( nb_realizations - 1 );
                variances[stat_id] = variance;
                // stddevs[stat_id] =std::sqrt( std::max( variance, 0.0 ) );
            }
        }

    public:
        std::vector< double > sum;
        std::vector< double > sum_squares;
        std::vector< double > means;
        std::vector< double > variances;
    };

} // namespace geode