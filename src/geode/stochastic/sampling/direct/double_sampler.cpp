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

#include <geode/stochastic/sampling/direct/double_sampler.hpp>

#include <geode/geometry/angle.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

namespace
{
    struct DistributionParams
    {
        double min;
        double max;
        double mean;
        double std;
    };

    DistributionParams complete_distribution_params(
        const geode::DoubleSampler::DistributionDescription& d )
    {
        DistributionParams p{};

        if( d.min_value && d.max_value )
        {
            p.min = *d.min_value;
            p.max = *d.max_value;
            p.mean = ( p.min + p.max ) / 2.0;
            p.std = ( p.max - p.min ) / std::sqrt( 12.0 );
            return p;
        }

        if( d.mean && d.standard_deviation )
        {
            p.mean = *d.mean;
            p.std = *d.standard_deviation;
            p.min = p.mean - std::sqrt( 3.0 ) * p.std;
            p.max = p.mean + std::sqrt( 3.0 ) * p.std;
            return p;
        }

        throw geode::OpenGeodeException(
            "[DistributionDescripption] Incomplete distribution description: "
            "need at least (min,max) or (mean,std)." );
        return p;
    }
} // namespace
namespace geode
{
    struct DistributionTypeHasher
    {
        std::size_t operator()( const DistributionType& d ) const noexcept
        {
            // Use the underlying string from NamedType
            return absl::Hash< std::string >{}( d.get() );
        }
    };

    using DistributionFactory = std::function< DoubleSampler::Distribution(
        const DoubleSampler::DistributionDescription& ) >;

    static absl::flat_hash_map< DistributionType, // key type
        DistributionFactory, // value type
        DistributionTypeHasher // custom hasher
        >
        distribution_registry = {
            { UniformClosed< double >::distribution_type_static(),
                []( const DoubleSampler::DistributionDescription& d ) {
                    auto p = complete_distribution_params( d );
                    UniformClosed< double > dist;
                    dist.min_value = p.min;
                    dist.max_value = p.max;
                    return dist;
                } },
            { UniformClosedOpen< double >::distribution_type_static(),
                []( const DoubleSampler::DistributionDescription& d ) {
                    auto p = complete_distribution_params( d );
                    UniformClosedOpen< double > dist;
                    dist.min_value = p.min;
                    dist.max_value = p.max;
                    return dist;
                } },
            { Gaussian::distribution_type_static(),
                []( const DoubleSampler::DistributionDescription& d ) {
                    auto p = complete_distribution_params( d );
                    Gaussian dist;
                    dist.mean = p.mean;
                    dist.standard_deviation = p.std;
                    return dist;
                } },
            { TruncatedGaussian::distribution_type_static(),
                []( const DoubleSampler::DistributionDescription& d ) {
                    auto p = complete_distribution_params( d );
                    TruncatedGaussian dist;
                    dist.mean = p.mean;
                    dist.standard_deviation = p.std;
                    dist.min_value = p.min;
                    dist.max_value = p.max;
                    return dist;
                } },
        };

    DoubleSampler::Distribution DoubleSampler::create_distribution(
        const DistributionDescription& desc )
    {
        auto it = distribution_registry.find( desc.distribution_type );
        if( it == distribution_registry.end() )
            throw geode::OpenGeodeException( absl::StrCat(
                "Unknown distribution type: ", desc.distribution_type.get() ) );
        return it->second( desc );
    }

    DoubleSampler::Distribution DoubleSampler::create_angle_distribution_in_rad(
        const DistributionDescription& desc_deg )
    {
        DistributionDescription desc_rad = desc_deg;
        if( desc_rad.mean )
        {
            auto mean_angle = Angle::create_from_degrees( *( desc_rad.mean ) );
            desc_rad.mean = mean_angle.normalized_between_0_and_2pi().radians();
        }
        if( desc_rad.standard_deviation )
        {
            auto std_angle =
                Angle::create_from_degrees( *( desc_rad.standard_deviation ) );
            desc_rad.standard_deviation =
                std_angle.normalized_between_0_and_2pi().radians();
        }
        if( desc_rad.min_value )
        {
            auto min_angle =
                Angle::create_from_degrees( *( desc_rad.min_value ) );
            desc_rad.min_value =
                min_angle.normalized_between_0_and_2pi().radians();
        }
        if( desc_rad.max_value )
        {
            auto max_angle =
                Angle::create_from_degrees( *( desc_rad.max_value ) );
            desc_rad.max_value =
                max_angle.normalized_between_0_and_2pi().radians();
        }
        // Call the general create_distribution
        return create_distribution( desc_rad );
    }

    double DoubleSampler::sample(
        RandomEngine& engine, const Distribution& dist )
    {
        return std::visit(
            [&engine]( auto&& d ) {
                using D = std::decay_t< decltype( d ) >;
                if constexpr( std::is_same_v< D, UniformClosed< double > > )
                    return engine.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, UniformClosedOpen< double > > )
                    return engine.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, Gaussian > )
                    return engine.sample_gaussian( d );
                if constexpr( std::is_same_v< D, TruncatedGaussian > )
                    return engine.sample_truncated_gaussian( d );
                throw OpenGeodeException( "DoubleSampler - Unsupported "
                                          "distribution for double" );
            },
            dist );
    }

} // namespace geode