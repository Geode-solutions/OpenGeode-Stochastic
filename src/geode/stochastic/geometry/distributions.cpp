/*
 * Copyright (c) 2019 - 2025 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <geode/stochastic/geometry/distributions.hpp>

#include <limits>

namespace geode
{
    template < typename Type >
    bool UniformClosed< Type >::is_valid() const
    {
        if( min_value < max_value )
        {
            return true;
        }
        if( min_value == max_value )
        {
            geode::Logger::warn(
                "[UniformClosed] - check boundaries definintion [", min_value,
                ",", max_value, "]." );
            return true;
        }
        geode::Logger::error(
            "[UniformClosed] - check boundaries definintion [", min_value, ",",
            max_value, "]." );
        return false;
    }
    template opengeode_stochastic_geometry_api struct UniformClosed< index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosed<
        local_index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosed<
        signed_index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosed< float >;
    template opengeode_stochastic_geometry_api struct UniformClosed< double >;

    template < typename Type >
    bool UniformClosedOpen< Type >::is_valid() const
    {
        if( min_value < max_value )
        {
            return true;
        }
        geode::Logger::error(
            "[UniformClosedOpen] - check boundaries definintion [", min_value,
            ",", max_value, "]." );
        return false;
    }
    template opengeode_stochastic_geometry_api struct UniformClosedOpen<
        index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosedOpen<
        local_index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosedOpen<
        signed_index_t >;
    template opengeode_stochastic_geometry_api struct UniformClosedOpen<
        float >;
    template opengeode_stochastic_geometry_api struct UniformClosedOpen<
        double >;

    bool Gaussian::is_valid() const
    {
        if( standard_deviation > 0 && std::isfinite( standard_deviation )
            && std::isfinite( mean ) )
        {
            return true;
        }
        geode::Logger::error(
            "[Gaussian] - check mean and standard deviation N(", mean, ",",
            standard_deviation, ")." );
        return false;
    }

    bool TruncatedGaussian::is_valid() const
    {
        if( standard_deviation <= 0 || std::isfinite( standard_deviation )
            || std::isfinite( mean ) )
        {
            geode::Logger::error(
                "[Gaussian] - check mean and standard deviation N(", mean, ",",
                standard_deviation, ")." );
            return false;
        }
        const auto max =
            max_value.value_or( std::numeric_limits< double >::infinity() );
        const auto min =
            min_value.value_or( -std::numeric_limits< double >::infinity() );

        if( min >= max )
        {
            geode::Logger::error( "[TruncatedGaussian] - check "
                                  "boundaries definintion [",
                min, ",", max, "]." );
            return false;
        }

        auto mean_out_of_bounds = ( mean < min ) || ( mean > max );
        if( mean_out_of_bounds )
        {
            geode::Logger::warn( "[TruncatedGaussian] - the mean (", mean,
                ") is defined out [", min, ",", max, "]." );
        }
        bool range_finite = std::isfinite( min ) && std::isfinite( max );
        bool small_range = false;
        if( range_finite )
        {
            small_range = ( std::abs( max - min ) < 0.9 * standard_deviation );
        }
        if( small_range )
        {
            geode::Logger::warn(
                "[TruncatedGaussian] - boundary range defined in [", min, ",",
                max, "] may be too narrow." );
        }
        return true;
    }
} // namespace geode
