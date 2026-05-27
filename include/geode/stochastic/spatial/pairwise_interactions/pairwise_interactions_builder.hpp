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

#include <geode/stochastic/common.hpp>

#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions_config.hpp>
#include <variant>

#pragma once
namespace geode
{
    template < typename ObjectType >
    std::unique_ptr< PairwiseInteraction< ObjectType > >
        build_pairwise_interaction_impl(
            const EuclideanDistanceCutoffConfig& cfg )
    {
        return std::make_unique< CenterEuclideanDistanceCutoff< ObjectType > >(
            cfg.threshold );
    }

    template < typename ObjectType >
    std::unique_ptr< PairwiseInteraction< ObjectType > >
        build_pairwise_interaction_impl(
            const MinimalDistanceCutoffConfig& cfg )
    {
        return std::make_unique< MinimalDistanceCutoff< ObjectType > >(
            cfg.threshold );
    }

    template < typename ObjectType, typename NewInteractionConfig >
    std::unique_ptr< PairwiseInteraction< ObjectType > >
        build_pairwise_interaction_impl(
            const NewInteractionConfig& /*unused*/ )
    {
        static_assert( sizeof( NewInteractionConfig ) == 0,
            "Unsupported PairwiseInteractionConfig type" );
        return nullptr;
    }

    template < typename ObjectType >
    std::unique_ptr< PairwiseInteraction< ObjectType > >
        build_pairwise_interaction_impl( const std::monostate& /*unused*/ )
    {
        throw OpenGeodeStochasticStochasticException{ nullptr,
            OpenGeodeException::TYPE::data,
            "[PairWiseInteractionBuilder] interaction config not initialized" };
    }

    template < typename ObjectType >
    std::unique_ptr< PairwiseInteraction< ObjectType > >
        build_pairwise_interaction( const PairwiseInteractionConfig& cfg )
    {
        return std::visit(
            []( auto&& interaction_cfg )
                -> std::unique_ptr< PairwiseInteraction< ObjectType > > {
                return build_pairwise_interaction_impl< ObjectType >(
                    interaction_cfg );
            },
            cfg );
    }
} // namespace geode