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

#include <geode/stochastic/applications/observed_intensity.hpp>

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>

namespace
{
    constexpr geode::index_t NB_MARK_SAMPLES{ 100'000 };
    constexpr std::string_view AREA_SEED{ "touching-domain-area" };

    // A fixed object outside the support of the length law of the set (all
    // the objects above a cutoff length are observed) is not a part of the
    // sampled population: it is added to it, as with free_objects.
    template < typename ObjectType >
    bool is_in_sampled_population( const ObjectType& object,
        const geode::ObjectSamplerConfig< ObjectType >& sampler_config )
    {
        if constexpr( !std::is_same_v< ObjectType, geode::OwnerSegment2D > )
        {
            geode_unused( object );
            geode_unused( sampler_config );
            return true;
        }
        else
        {
            const auto& length_law = sampler_config.length;
            const auto length = object.length();
            if( length_law.max_value
                && length
                       > length_law.max_value.value() + geode::GLOBAL_EPSILON )
            {
                return false;
            }
            return !length_law.min_value
                   || length >= length_law.min_value.value()
                                    - geode::GLOBAL_EPSILON;
        }
    }
} // namespace

namespace geode
{
    template < typename ObjectType >
    double touching_domain_area( const SpatialDomain< ObjectType::dim >& domain,
        const ObjectSamplerConfig< ObjectType >& sampler_config )
    {
        if constexpr( !std::is_same_v< ObjectType, OwnerSegment2D > )
        {
            geode_unused( sampler_config );
            return domain.n_volume();
        }
        else
        {
            // For a segment of length L and azimuth theta (direction
            // (sin theta, cos theta)), the anchors for which it touches the
            // w x h domain form the Minkowski sum of the domain and of the
            // reversed segment: w * h + L * (w |cos theta| + h |sin theta|).
            // The buffer, at least the maximal length, contains it. Only the
            // marks are sampled: the precision does not depend on the size of
            // the domain.
            const auto length =
                DoubleSampler::create_distribution( sampler_config.length );
            const auto azimuth =
                DoubleSampler::create_distribution( sampler_config.azimuth );
            const auto& box = domain.box();
            const auto width = box.max().value( 0 ) - box.min().value( 0 );
            const auto height = box.max().value( 1 ) - box.min().value( 1 );
            RandomEngine engine;
            engine.set_seed( AREA_SEED );
            double boundary_area{ 0. };
            for( const auto sample_id : Range{ NB_MARK_SAMPLES } )
            {
                geode_unused( sample_id );
                const auto segment_length =
                    DoubleSampler::sample( engine, length );
                const auto segment_azimuth =
                    DoubleSampler::sample( engine, azimuth );
                boundary_area +=
                    segment_length
                    * ( width * std::abs( std::cos( segment_azimuth ) )
                        + height * std::abs( std::sin( segment_azimuth ) ) );
            }
            return width * height + boundary_area / NB_MARK_SAMPLES;
        }
    }

    template < typename ObjectType >
    void add_observation_term( SimulationContextConfig< ObjectType >& config,
        index_t set_index,
        double lambda )
    {
        const auto& set = config.sets[set_index];
        const auto domain = build_spatial_domain( config.domain );
        double nb_fixed{ 0 };
        index_t nb_outside_law{ 0 };
        for( const auto& object : set.fixed_objects )
        {
            if( !is_in_sampled_population( object, set.sampler ) )
            {
                nb_outside_law++;
                continue;
            }
            if( SpatialDomainChecker< ObjectType >::touches_domain(
                    *domain, object ) )
            {
                nb_fixed++;
            }
        }
        if( nb_outside_law > 0 )
        {
            Logger::info( "[IntensityScope] Object set '", set.name,
                "': ", nb_outside_law,
                " fixed objects outside the length law are added to the "
                "simulated ones, not counted as observed ones" );
        }
        const auto nb_expected =
            lambda * touching_domain_area( *domain, set.sampler );
        OpenGeodeStochasticStochasticException::check_exception(
            nb_expected - nb_fixed >= -GLOBAL_EPSILON, nullptr,
            OpenGeodeException::TYPE::data, "[IntensityScope] Object set '",
            set.name, "': lambda = ", lambda, " gives ", nb_expected,
            " objects expected to touch the domain, fewer than its ", nb_fixed,
            " fixed objects touching it: lambda is underestimated." );
        const auto observed_ratio =
            nb_expected > 0. ? std::min( nb_fixed / nb_expected, 1. ) : 0.;
        Logger::info( "[IntensityScope] Object set '", set.name,
            "': ", nb_fixed, " fixed objects touch the domain, out of ",
            nb_expected,
            " expected: free objects touching the domain have the density ",
            lambda * ( 1. - observed_ratio ) );

        SingleObjectTermConfig observation;
        observation.term_name = config.model.unique_term_name(
            absl::StrCat( set.name, "_observation" ) );
        observation.object_set_names = { set.name };
        observation.lambda = 1. - observed_ratio;
        observation.object_feature = ObjectTouchingDomainFeatureConfig{};
        config.model.terms.emplace_back( std::move( observation ) );
    }

    template opengeode_stochastic_stochastic_api double touching_domain_area(
        const SpatialDomain< 2 >&, const ObjectSamplerConfig< Point2D >& );
    template opengeode_stochastic_stochastic_api double touching_domain_area(
        const SpatialDomain< 3 >&, const ObjectSamplerConfig< Point3D >& );
    template opengeode_stochastic_stochastic_api double touching_domain_area(
        const SpatialDomain< 2 >&,
        const ObjectSamplerConfig< OwnerSegment2D >& );

    template opengeode_stochastic_stochastic_api void add_observation_term(
        SimulationContextConfig< Point2D >&, index_t, double );
    template opengeode_stochastic_stochastic_api void add_observation_term(
        SimulationContextConfig< Point3D >&, index_t, double );
    template opengeode_stochastic_stochastic_api void add_observation_term(
        SimulationContextConfig< OwnerSegment2D >&, index_t, double );
} // namespace geode
