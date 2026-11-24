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

#include <geode/stochastic/common.hpp>

#include <cmath>
#include <limits>
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
            // Equal bounds are valid: every sample is this value.
            if( std::isfinite( min_value ) && std::isfinite( max_value )
                && min_value <= max_value )
            {
                return true;
            }
            geode::Logger::error(
                "[Uniform Closed] - check range boundaries definition [",
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
            return absl::StrCat( distribution_type().get(), "[", min_value, ",",
                max_value, "]" );
        }
    };

    template < typename Type >
    struct UniformClosedOpen
    {
        UniformClosedOpen() = default;
        bool is_valid() const
        {
            if( std::isfinite( min_value ) && std::isfinite( max_value )
                && min_value < max_value )
            {
                return true;
            }
            geode::Logger::error(
                "[Uniform ClosedOpen] - check range boundaries definition [",
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
            return absl::StrCat( distribution_type().get(), "[", min_value, ",",
                max_value, ")" );
        }
    };

    struct opengeode_stochastic_stochastic_api Gaussian
    {
        Gaussian() = default;
        bool is_valid() const
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

        bool is_valid() const
        {
            if( standard_deviation <= 0 || !std::isfinite( standard_deviation )
                || !std::isfinite( mean ) )
            {
                geode::Logger::error( "[Truncated Gaussian] - check mean and "
                                      "standard deviation N(",
                    mean, ",", standard_deviation, ")." );
                return false;
            }
            const auto max =
                max_value.value_or( std::numeric_limits< double >::infinity() );
            const auto min = min_value.value_or(
                -std::numeric_limits< double >::infinity() );

            // written in positive form to reject NaN bounds
            if( !( min < max ) )
            {
                geode::Logger::error( "[Truncated Gaussian] - check "
                                      "range boundaries definition [",
                    min, ",", max, "]." );
                return false;
            }

            if( min > mean + 6.0 * standard_deviation
                || max < mean - 6.0 * standard_deviation )
            {
                geode::Logger::error( "[Truncated Gaussian] - Truncation "
                                      "range is too far from mean. N(",
                    mean, ",", standard_deviation, ") is defined in [", min,
                    ",", max, "]." );
                return false;
            }
            return true;
        }

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
            return absl::StrCat( distribution_type().get(), "(", mean, ",",
                standard_deviation, ") in [", min_str, ",", max_str, "]" );
        }
    };

    struct opengeode_stochastic_stochastic_api VonMises
    {
        VonMises() = default;

        // 2D von Mises parameters
        double mean{ 0.0 }; // mean direction in radians
        double concentration{ 1.0 }; // kappa (concentration parameter)

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "VonMises" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        [[nodiscard]] bool is_valid() const
        {
            if( concentration >= 0.0 && std::isfinite( concentration )
                && std::isfinite( mean ) )
            {
                return true;
            }
            geode::Logger::error( "[VonMises] - check mean and concentration (",
                mean, ",", concentration, ")." );
            return false;
        }

        std::string string() const
        {
            return absl::StrCat( distribution_type().get(), "(mean=", mean,
                ", kappa=", concentration, ")" );
        }
    };

    struct TruncatedLogNormal
    {
        TruncatedLogNormal() = default;

        // Parameters of the underlying normal distribution
        double mean{ 0.0 };
        double standard_deviation{ 1.0 };

        std::optional< double > min_value;
        std::optional< double > max_value;

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "TruncatedLogNormal" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        bool is_valid() const
        {
            if( standard_deviation <= 0 || !std::isfinite( standard_deviation )
                || !std::isfinite( mean ) )
            {
                geode::Logger::error( "[TruncatedLogNormal] - check mean and "
                                      "standard deviation N(",
                    mean, ",", standard_deviation, ")." );
                return false;
            }
            const auto max =
                max_value.value_or( std::numeric_limits< double >::infinity() );
            const auto min = min_value.value_or( 0. );

            if( min < 0. )
            {
                geode::Logger::error( "[TruncatedLogNormal] - the "
                                      "lower bound cannot be negative: ",
                    min, "." );
                return false;
            }
            if( !( min < max ) )
            {
                geode::Logger::error( "[TruncatedLogNormal] - check "
                                      "range boundaries definition [",
                    min, ",", max, "]." );
                return false;
            }
            return true;
        }

        std::string string() const
        {
            std::string min_str = min_value.has_value()
                                      ? std::to_string( min_value.value() )
                                      : "0.";

            std::string max_str = max_value.has_value()
                                      ? std::to_string( max_value.value() )
                                      : "+inf";
            return absl::StrCat( distribution_type().get(), "(", mean, ",",
                standard_deviation, ") in [", min_str, ",", max_str, "]" );
        }
    };

    struct TruncatedPowerLaw
    {
        TruncatedPowerLaw() = default;

        // Power-law exponent
        double alpha{ 2.0 }; // default value > 0

        // Optional truncation bounds
        std::optional< double > min_value;
        std::optional< double > max_value;

        [[nodiscard]] static DistributionType distribution_type_static()
        {
            return DistributionType{ "TruncatedPowerLaw" };
        }

        [[nodiscard]] DistributionType distribution_type() const
        {
            return distribution_type_static();
        }

        // Lower bound used when min_value is not given
        [[nodiscard]] double lower_bound() const
        {
            return min_value.value_or( GLOBAL_EPSILON );
        }

        // Check if the parameters are valid
        bool is_valid() const
        {
            if( alpha <= 0.0 || !std::isfinite( alpha ) )
            {
                geode::Logger::error(
                    "[TruncatedPowerLaw] - exponent alpha must be > 0: alpha=",
                    alpha );
                return false;
            }

            const double xmin = lower_bound();
            const double xmax =
                max_value.value_or( std::numeric_limits< double >::infinity() );

            if( xmin <= 0.0 )
            {
                geode::Logger::error(
                    "[TruncatedPowerLaw] - min_value must be positive: ",
                    xmin );
                return false;
            }

            if( !( xmin < xmax ) )
            {
                geode::Logger::error(
                    "[TruncatedPowerLaw] - check range boundaries definition [",
                    xmin, ",", xmax, "]" );
                return false;
            }

            if( alpha <= 1.0 && !std::isfinite( xmax ) )
            {
                geode::Logger::error( "[TruncatedPowerLaw] - with alpha <= 1 "
                                      "(here ",
                    alpha, "), the law needs a finite max_value." );
                return false;
            }

            return true;
        }

        // String representation
        std::string string() const
        {
            const std::string min_str =
                min_value.has_value() ? std::to_string( min_value.value() )
                                      : "ε";
            const std::string max_str =
                max_value.has_value() ? std::to_string( max_value.value() )
                                      : "+inf";

            return absl::StrCat( distribution_type().get(), "(alpha=", alpha,
                ") in [", min_str, ",", max_str, "]" );
        }
    };
} // namespace geode