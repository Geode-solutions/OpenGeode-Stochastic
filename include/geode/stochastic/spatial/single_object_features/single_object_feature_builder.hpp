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

#include <variant>

#include <geode/basic/assert.hpp>

#include <geode/stochastic/spatial/single_object_features/segment_length_feature.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>

namespace geode
{
    template < typename ObjectType >
    std::unique_ptr< SingleObjectFeature< ObjectType > >
        build_single_object_feature_impl(
            const ObjectInDomainFeatureConfig& cfg )
    {
        geode_unused( cfg );
        return std::make_unique< ObjectInDomainFeature< ObjectType > >();
    };

    template < typename ObjectType >
    std::unique_ptr< SingleObjectFeature< ObjectType > >
        build_single_object_feature_impl(
            const SegmentLengthInsideBoxFeatureConfig& cfg )
    {
        if constexpr( std::is_same_v< ObjectType, OwnerSegment2D > )
        {
            return std::make_unique< SegmentLengthInsideBoxFeature >(
                cfg.characteristic_length );
        }
        else
        {
            throw OpenGeodeStochasticStochasticException{ nullptr,
                OpenGeodeException::TYPE::data,
                "[SingleObjectFeatureBuilder] SegmentLengthInsideBoxFeature "
                "not valid for this ObjectType" };
        }
    }

    template < typename ObjectType, typename NewObjectFeatureConfig >
    std::unique_ptr< SingleObjectFeature< ObjectType > >
        build_single_object_feature_impl(
            const NewObjectFeatureConfig& /*unused*/ )
    {
        static_assert( sizeof( NewObjectFeatureConfig ) == 0,
            "[SingleObjectFeatureBuilder] Unsupported "
            "SingleObjectFeatureConfig type" );
        return nullptr;
    }

    template < typename ObjectType >
    std::unique_ptr< SingleObjectFeature< ObjectType > >
        build_single_object_feature_impl( const std::monostate& /*unused*/ )
    {
        throw OpenGeodeStochasticStochasticException{ nullptr,
            OpenGeodeException::TYPE::data,
            "[SingleObjectFeatureBuilder] object feature "
            "config not initialized" };
    }

    template < typename ObjectType >
    std::unique_ptr< SingleObjectFeature< ObjectType > >
        build_single_object_feature( const SingleObjectFeatureConfig& cfg )
    {
        return std::visit(
            []( auto&& interaction_cfg )
                -> std::unique_ptr< SingleObjectFeature< ObjectType > > {
                return build_single_object_feature_impl< ObjectType >(
                    interaction_cfg );
            },
            cfg );
    }

} // namespace geode