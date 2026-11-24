
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

#include <geode/stochastic/sampling/distributions.hpp>

#include <variant>

namespace geode
{
    class RandomEngine;
} // namespace geode

namespace geode
{

    struct opengeode_stochastic_stochastic_api DoubleSampler
    {
        using Distribution = std::variant< UniformClosed< double >,
            UniformClosedOpen< double >,
            Gaussian,
            TruncatedGaussian,
            VonMises,
            TruncatedLogNormal,
            TruncatedPowerLaw >;

        struct DistributionDescription
        {
            explicit DistributionDescription( absl::string_view name_in )
                : name{ name_in }
            {
            }

            DistributionDescription() = default;

            std::string name{ "default_distribution" };
            DistributionType distribution_type;

            std::optional< double > min_value;
            std::optional< double > max_value;
            std::optional< double > mean;
            std::optional< double > standard_deviation;
            std::optional< double >
                kappa; // concentration parameter for VonMises distribution
            std::optional< double >
                alpha; // exponent parameter for power law distribution

            [[nodiscard]] std::string string() const
            {
                auto message = absl::StrCat( "Distribution - ", name );
                absl::StrAppend( &message,
                    "\n\t - distribution type: ", distribution_type.get() );
                if( min_value.has_value() )
                {
                    absl::StrAppend(
                        &message, "\n\t - min value: ", min_value.value() );
                }
                if( max_value.has_value() )
                {
                    absl::StrAppend(
                        &message, "\n\t - max value: ", max_value.value() );
                }
                if( mean.has_value() )
                {
                    absl::StrAppend(
                        &message, "\n\t - mean value: ", mean.value() );
                }
                if( standard_deviation.has_value() )
                {
                    absl::StrAppend( &message,
                        "\n\t - std value: ", standard_deviation.value() );
                }
                if( kappa.has_value() )
                {
                    absl::StrAppend( &message,
                        "\n\t - kappa (von mises concentration) value: ",
                        kappa.value() );
                }
                if( alpha.has_value() )
                {
                    absl::StrAppend( &message,
                        "\n\t - alpha (power law exponent) value: ",
                        alpha.value() );
                }
                return message;
            }
        };

        static Distribution create_distribution(
            const DistributionDescription& desc );

        /*!
         * Creates an angle distribution in radians from a description in
         * degrees. Bounds, mean and standard deviation are converted without
         * normalization: [-90, 90] gives [-pi/2, pi/2] and [0, 360] the full
         * circle. Only the mean direction of a von Mises law is normalized in
         * [0, 2pi[. The log-normal mean and standard deviation describe
         * ln(x): the mean is shifted by ln(pi/180), the deviation is kept.
         * Use UniformClosedOpen on [0, 360] to cover the circle
         * without drawing both 0 and 2pi. Samples are not normalized: use
         * Angle::normalized_between_0_and_2pi() on them if needed.
         */
        static Distribution create_rad_angle_distribution_from_degree(
            const DistributionDescription& angle_desc_deg );

        static double sample( RandomEngine& engine, const Distribution& dist );
    };

} // namespace geode