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

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

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
        const DistributionDescription& ) >;

    static absl::flat_hash_map< DistributionType, // key type
        DistributionFactory, // value type
        DistributionTypeHasher // custom hasher
        >
        distribution_registry = {
            { UniformClosed< double >::distribution_type_static(),
                []( const DistributionDescription& d ) {
                    UniformClosed< double > dist;
                    dist.min_value = d.min_value;
                    dist.max_value = d.max_value;
                    return dist;
                } },
            { UniformClosedOpen< double >::distribution_type_static(),
                []( const DistributionDescription& d ) {
                    UniformClosedOpen< double > dist;
                    dist.min_value = d.min_value;
                    dist.max_value = d.max_value;
                    return dist;
                } },
            { Gaussian::distribution_type_static(),
                []( const DistributionDescription& d ) {
                    Gaussian dist;
                    dist.mean =
                        d.mean.value_or( ( d.max_value - d.min_value ) / 2.0 );
                    dist.standard_deviation = d.standard_deviation.value_or(
                        ( d.max_value - d.min_value ) / 6.0 );
                    return dist;
                } },
            { TruncatedGaussian::distribution_type_static(),
                []( const DistributionDescription& d ) {
                    TruncatedGaussian dist;
                    dist.mean =
                        d.mean.value_or( ( d.max_value - d.min_value ) / 2.0 );
                    dist.standard_deviation = d.standard_deviation.value_or(
                        ( d.max_value - d.min_value ) / 6.0 );
                    dist.min_value = d.min_value;
                    dist.max_value = d.max_value;
                    return dist;
                } }
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