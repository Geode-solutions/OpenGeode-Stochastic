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

#include <variant>

#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>

#pragma once
namespace geode
{
    struct EuclideanDistanceCutoffConfig
    {
        double threshold;
    };

    struct MinimalDistanceCutoffConfig
    {
        double threshold;
    };

    using PairwiseInteractionConfig =
        std::variant< EuclideanDistanceCutoffConfig,
            MinimalDistanceCutoffConfig >;

    template < typename ObjectType >
    std::unique_ptr< PairwiseInteraction< ObjectType > > build_interaction(
        const PairwiseInteractionConfig& cfg )
    {
        return std::visit(
            [&]( auto&& c )
                -> std::unique_ptr< PairwiseInteraction< ObjectType > > {
                using T = std::decay_t< decltype( c ) >;

                if constexpr( std::is_same_v< T,
                                  EuclideanDistanceCutoffConfig > )
                {
                    return std::make_unique<
                        CenterEuclideanDistanceCutoff< ObjectType > >(
                        c.threshold );
                }
                else if constexpr( std::is_same_v< T,
                                       MinimalDistanceCutoffConfig > )
                {
                    return std::make_unique<
                        MinimalDistanceCutoff< ObjectType > >(
                        c.threshold, c.weight );
                }
            },
            cfg );
    }
} // namespace geode