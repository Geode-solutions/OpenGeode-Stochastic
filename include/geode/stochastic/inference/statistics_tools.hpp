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
#include <geode/stochastic/inference/target_statistics.hpp>

namespace geode
{
    namespace statistics
    {
        template < typename ObjectType >
        void validate( const StatisticsTracker< ObjectType >& tracker,
            const TargetStatistics< ObjectType >& targets )
        {
            const auto& model = targets.model();

            for( const auto& term_uuid : targets.active_terms() )
            {
                const auto mean = tracker.mean( term_uuid );
                const auto target = targets.target( term_uuid );

                const auto rel_error =
                    std::fabs( mean - target )
                    / ( std::fabs( target ) + geode::GLOBAL_EPSILON );

                OpenGeodeStochasticStochasticException::check_exception(
                    rel_error < targets.tolerance( term_uuid ), nullptr,
                    OpenGeodeException::TYPE::result,
                    "[StatisticsValidator] Failure for term ",
                    model.term_name( term_uuid ), "\n  mean    = ", mean,
                    "\n  target  = ", target, "\n  error   = ", rel_error,
                    "\n  tol     = ", targets.tolerance( term_uuid ) );
            }
        }

        template < typename ObjectType >
        double quadratic_loss( const StatisticsTracker< ObjectType >& tracker,
            const TargetStatistics< ObjectType >& targets )
        {
            double loss = 0.0;

            for( const auto& term_uuid : targets.active_terms() )
            {
                const auto diff =
                    tracker.mean( term_uuid ) - targets.value( term_uuid );

                loss += diff * diff;
            }

            return loss;
        }
    } // namespace statistics
} // namespace geode