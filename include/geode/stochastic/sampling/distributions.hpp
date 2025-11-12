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

#include <optional>

#include <geode/basic/named_type.hpp>

namespace geode
{
    struct DistributionTag
    {
    };

    using DistributionType = NamedType< std::string, DistributionTag >;

    template < typename Type >
    struct UniformClosed
    {
        UniformClosed() = default;
        bool is_valid() const
        {
            if( min_value < max_value )
            {
                return true;
            }
            if( min_value == max_value )
            {
                geode::Logger::warn(
                    "[Uniform Closed] - check range boundaries definintion [",
                    min_value, ",", max_value, "]." );
                return true;
            }
            geode::Logger::error(
                "[Uniform Closed] - check range boundaries definintion [",
                min_value, ",", max_value, "]." );
            return false;
        }

        Type min_value{ static_cast< Type >( 0 ) };
        Type max_value{ static_cast< Type >( 1 ) };

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "UniformClosed" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }
        std::string string() const
        {
            return absl::StrCat(
                distribution_type().get(), "[", min_value, max_value, "]" );
        }
    };

    template < typename Type >
    struct UniformClosedOpen
    {
        UniformClosedOpen() = default;
        bool is_valid() const
        {
            if( min_value < max_value )
            {
                return true;
            }
            geode::Logger::error(
                "[Uniform ClosedOpen] - check range boundaries definintion [",
                min_value, ",", max_value, "]." );
            return false;
        }

        Type min_value{ static_cast< Type >( 0 ) };
        Type max_value{ static_cast< Type >( 1 ) };

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "UniformClosedOpen" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        std::string string() const
        {
            return absl::StrCat(
                distribution_type().get(), "[", min_value, max_value, "]" );
        }
    };

    struct opengeode_stochastic_stochastic_api Gaussian
    {
        Gaussian() = default;
        bool is_valid() const;

        double mean{ 0. };
        double standard_deviation{ 1. };

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "Gaussian" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        std::string string() const
        {
            return absl::StrCat( distribution_type().get(), "(", mean, ",",
                standard_deviation, ")" );
        }
    };

    struct opengeode_stochastic_stochastic_api TruncatedGaussian
    {
        TruncatedGaussian() = default;

        bool is_valid() const;

        double mean{ 0. };
        double standard_deviation{ 1. };

        std::optional< double > min_value;
        std::optional< double > max_value;

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "TruncatedGaussian" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        std::string string() const
        {
            std::string min_str = min_value.has_value()
                                      ? std::to_string( min_value.value() )
                                      : "-inf";

            std::string max_str = max_value.has_value()
                                      ? std::to_string( max_value.value() )
                                      : "+inf";
            return absl::StrCat( distribution_type().get(), "(", mean,
                standard_deviation, ") in [", min_str, ",", max_str, "]" );
        }
    };

} // namespace geode