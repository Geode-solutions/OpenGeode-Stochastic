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

#include <geode/stochastic/geometry/sampler/double_distribution_sampler.hpp>

#include <geode/stochastic/geometry/common.hpp>
#include <geode/stochastic/geometry/random_engine.hpp>

namespace geode
{
    double DistributionSampler::sample( const DoubleDistribution& dist )
    {
        return std::visit(
            [this]( auto&& d ) {
                using D = std::decay_t< decltype( d ) >;
                if constexpr( std::is_same_v< D, UniformClosed< double > > )
                    return engine_.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, UniformClosedOpen< double > > )
                    return engine_.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, Gaussian > )
                    return engine_.sample_gaussian( d );
                if constexpr( std::is_same_v< D, TruncatedGaussian > )
                    return engine_.sample_truncated_gaussian( d );
                OPENGEODE_EXCEPTION(
                    "RandomSampler - Unsupported distribution for double" );
            },
            dist );
    }
} // namespace geode