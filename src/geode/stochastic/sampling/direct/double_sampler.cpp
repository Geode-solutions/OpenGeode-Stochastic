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
#include <cmath>

#include <absl/container/flat_hash_map.h>

#include <geode/stochastic/sampling/direct/double_sampler.hpp>

#include <geode/geometry/angle.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

namespace geode
{
    struct DistributionTypeHasher
    {
        std::size_t operator()( const DistributionType& distrib ) const noexcept
        {
            // Use the underlying string from NamedType
            return absl::Hash< std::string >{}( distrib.get() );
        }
    };

    using DistributionFactory = std::function< DoubleSampler::Distribution(
        const DoubleSampler::DistributionDescription& ) >;

    static const auto& distribution_registry()
    {
        static const absl::flat_hash_map< DistributionType, // key type
            DistributionFactory, // value type
            DistributionTypeHasher // custom hasher
            >
            registry = {
                { UniformClosed< double >::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.min_value && desc.max_value, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "Uniform distribution need at least min "
                            "and max values" );
                        UniformClosed< double > dist;
                        dist.min_value = desc.min_value.value();
                        dist.max_value = desc.max_value.value();
                        if( dist.min_value == dist.max_value )
                        {
                            Logger::warning( "[DoubleSampler] Uniform "
                                             "distribution with equal bounds: "
                                             "every sample is ",
                                dist.min_value, "." );
                        }
                        return dist;
                    } },
                { UniformClosedOpen< double >::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.min_value && desc.max_value, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "Uniform distribution need at least min "
                            "and max values" );
                        UniformClosedOpen< double > dist;
                        dist.min_value = desc.min_value.value();
                        dist.max_value = desc.max_value.value();
                        return dist;
                    } },
                { Gaussian::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.mean && desc.standard_deviation, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "Gaussian distribution need at least mean "
                            "and standard deviation values" );
                        Gaussian dist;
                        dist.mean = desc.mean.value();
                        dist.standard_deviation =
                            desc.standard_deviation.value();
                        return dist;
                    } },
                { TruncatedGaussian::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.mean && desc.standard_deviation, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "Truncated Gaussian distribution need at least "
                            "mean "
                            "and standard deviation values" );
                        TruncatedGaussian dist;
                        dist.mean = desc.mean.value();
                        dist.standard_deviation =
                            desc.standard_deviation.value();
                        dist.min_value = desc.min_value;
                        dist.max_value = desc.max_value;

                        return dist;
                    } },
                { VonMises::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.mean && desc.kappa, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "Von Mises distribution need at least mean "
                            "and concentration (kappa) values" );
                        VonMises dist;
                        dist.mean = desc.mean.value();
                        dist.concentration = desc.kappa.value();
                        return dist;
                    } },
                { TruncatedLogNormal::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.mean && desc.standard_deviation, nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "TruncatedLogNormal distribution need mean "
                            "and standard deviation values of the underlying "
                            "normal distribution." );
                        TruncatedLogNormal dist;
                        dist.mean = desc.mean.value();
                        dist.standard_deviation =
                            desc.standard_deviation.value();
                        dist.min_value = desc.min_value;
                        dist.max_value = desc.max_value;
                        return dist;
                    } },
                { TruncatedPowerLaw::distribution_type_static(),
                    []( const DoubleSampler::DistributionDescription& desc ) {
                        OpenGeodeStochasticStochasticException::check_exception(
                            desc.alpha.has_value(), nullptr,
                            OpenGeodeException::TYPE::data,
                            "[DoubleSampler] Incomplete description for "
                            "TruncatedPowerLaw distribution need power law "
                            "exponent (alpha)." );
                        TruncatedPowerLaw dist;
                        dist.alpha = desc.alpha.value();
                        dist.min_value = desc.min_value;
                        dist.max_value = desc.max_value;
                        return dist;
                    } },
            };
        return registry;
    }

    DoubleSampler::Distribution DoubleSampler::create_distribution(
        const DistributionDescription& desc )
    {
        auto it = distribution_registry().find( desc.distribution_type );
        if( it == distribution_registry().end() )
        {
            throw OpenGeodeStochasticStochasticException{ nullptr,
                OpenGeodeException::TYPE::data,
                "[DoubleSampler] Unknown distribution type: ",
                desc.distribution_type.get() };
        }
        auto distribution = it->second( desc );
        std::visit(
            [&desc]( const auto& law ) {
                OpenGeodeStochasticStochasticException::check_exception(
                    law.is_valid(), nullptr, OpenGeodeException::TYPE::data,
                    "[DoubleSampler] Invalid distribution: ", desc.string() );
            },
            distribution );
        return distribution;
    }

    DoubleSampler::Distribution
        DoubleSampler::create_rad_angle_distribution_from_degree(
            const DistributionDescription& desc_deg )
    {
        // Degrees are converted to radians without normalization: 0-360
        // stays the full circle, -90-90 a range across 0, and a standard
        // deviation keeps its value. Only the mean direction of a von Mises
        // law is normalized in [0, 2pi[.
        const auto to_radians = []( std::optional< double >& value ) {
            if( value )
            {
                value = Angle::create_from_degrees( value.value() ).radians();
            }
        };
        DistributionDescription desc_rad = desc_deg;
        to_radians( desc_rad.min_value );
        to_radians( desc_rad.max_value );
        if( desc_rad.distribution_type
            == TruncatedLogNormal::distribution_type_static() )
        {
            // The parameters describe ln(x): scaling x by pi/180 shifts the
            // mean by ln(pi/180) and keeps the standard deviation.
            if( desc_rad.mean )
            {
                desc_rad.mean.value() += std::log( M_PI / 180. );
            }
            return create_distribution( desc_rad );
        }
        to_radians( desc_rad.mean );
        to_radians( desc_rad.standard_deviation );
        if( desc_rad.mean
            && desc_rad.distribution_type
                   == VonMises::distribution_type_static() )
        {
            desc_rad.mean = Angle::create_from_radians( desc_rad.mean.value() )
                                .normalized_between_0_and_2pi()
                                .radians();
        }
        return create_distribution( desc_rad );
    }

    double DoubleSampler::sample(
        RandomEngine& engine, const Distribution& dist )
    {
        return std::visit(
            [&engine]( auto&& cur_dist ) {
                using D = std::decay_t< decltype( cur_dist ) >;
                if constexpr( std::is_same_v< D, UniformClosed< double > > )
                {
                    return engine.sample_uniform< double >( cur_dist );
                }
                if constexpr( std::is_same_v< D, UniformClosedOpen< double > > )
                {
                    return engine.sample_uniform< double >( cur_dist );
                }
                if constexpr( std::is_same_v< D, Gaussian > )
                {
                    return engine.sample_gaussian( cur_dist );
                }
                if constexpr( std::is_same_v< D, TruncatedGaussian > )
                {
                    return engine.sample_truncated_gaussian( cur_dist );
                }
                if constexpr( std::is_same_v< D, VonMises > )
                {
                    return engine.sample_von_mises( cur_dist );
                }
                if constexpr( std::is_same_v< D, TruncatedLogNormal > )
                {
                    return engine.sample_truncated_lognormal( cur_dist );
                }
                if constexpr( std::is_same_v< D, TruncatedPowerLaw > )
                {
                    return engine.sample_truncated_powerlaw( cur_dist );
                }
                throw OpenGeodeStochasticStochasticException{ nullptr,
                    OpenGeodeException::TYPE::data,
                    "DoubleSampler - Unsupported distribution for double" };
            },
            dist );
    }

} // namespace geode