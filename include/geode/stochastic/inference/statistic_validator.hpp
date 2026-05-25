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

#include <absl/container/flat_hash_map.h>

#include <geode/stochastic/inference/statistics_tracker.hpp>
#include <geode/stochastic/inference/target_statistic.hpp>

namespace geode
{
    template < typename ObjectType >
    class StatisticsValidator
    {
    public:
        void check( const StatisticsTracker& monitor,
            const std::vector< TargetStatistic >& targets ) const
        {
            for( const auto& target : targets )
            {
                const auto mean = monitor.mean( target.term_id );
                const auto rel_error = std::fabs( mean - target.value )
                                       / ( std::fabs( target.value ) + 1e-12 );
                OpenGeodeStochasticStochasticException::check_exception(
                    rel_error < target.tolerance, nullptr,
                    OpenGeodeException::TYPE::result,
                    "[StatisticsValidator] Failure for term ",
                    t.term_id.string(), "\n  mean    = ", mean,
                    "\n  target  = ", target.value,
                    "\n  error
                    = ", rel_error,
                      "\n  tol     = ",
                    target.tolerance );
            }
        }
    };
} // namespace geode