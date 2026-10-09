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

#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>

#include <geode/geometry/angle.hpp>
#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/basic_objects/sphere.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>
#include <geode/stochastic/sampling/direct/segment_uniform_sampler.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace
{
    // Segments anchored outside the extended domain are not simulated: the
    // buffer must be larger than the longest segment so that every segment
    // crossing the domain can be generated (no edge bias on the statistics).
    void check_buffer_size( const geode::SpatialDomain< 2 >& domain,
        const geode::DoubleSampler::DistributionDescription& length )
    {
        const auto bounded =
            length.max_value.has_value()
            && length.distribution_type
                   != geode::Gaussian::distribution_type_static();
        if( !bounded )
        {
            geode::Logger::warn( "[UniformSegmentSetSampler] Unbounded "
                                 "length distribution: segments longer than "
                                 "the buffer size (",
                domain.buffer_size(),
                ") bias the statistics near the domain boundary." );
            return;
        }
        geode::OpenGeodeStochasticStochasticException::check_exception(
            length.max_value.value() <= domain.buffer_size(), nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[UniformSegmentSetSampler] Buffer size (", domain.buffer_size(),
            ") must be larger than the maximal segment length (",
            length.max_value.value(),
            "): otherwise segments anchored outside the extended domain and "
            "crossing the domain are missing." );
    }

    geode::Vector2D direction_from_azimuth( double azimuth_rad )
    {
        const auto azimuth = geode::Angle::create_from_radians( azimuth_rad );
        return geode::Vector2D( { azimuth.sin(), azimuth.cos() } );
    }
} // namespace

namespace geode
{
    UniformSegmentSetSampler::UniformSegmentSetSampler(
        const SpatialDomain< 2 >& domain,
        const ObjectSamplerConfig< OwnerSegment2D >& config )
        : domain_{ domain },
          length_{ DoubleSampler::create_distribution( config.length ) },
          azimuth_{ DoubleSampler::create_distribution( config.azimuth ) },
          move_ratio_{ config.move_ratio }
    {
        check_buffer_size( domain_, config.length );
        OpenGeodeStochasticStochasticException::check_exception(
            config.translation_ratio >= 0. && config.rotation_ratio >= 0.
                && config.stretch_ratio >= 0.,
            nullptr, OpenGeodeException::TYPE::data,
            "[UniformSegmentSetSampler] Change move ratios cannot be "
            "negative." );
        const auto total_ratio = config.translation_ratio
                                 + config.rotation_ratio + config.stretch_ratio;
        OpenGeodeStochasticStochasticException::check_exception(
            total_ratio > 0., nullptr, OpenGeodeException::TYPE::data,
            "[UniformSegmentSetSampler] At least one change move ratio must be "
            "positive." );
        translation_probability_ = config.translation_ratio / total_ratio;
        const auto mark_ratio = config.rotation_ratio + config.stretch_ratio;
        // probability of a rotation knowing that the move is not a translation
        rotation_probability_ =
            mark_ratio > 0. ? config.rotation_ratio / mark_ratio : 0.;

        auto volume = domain_.extended_n_volume();
        OpenGeodeStochasticStochasticException::check_exception( volume != 0.,
            nullptr, OpenGeodeException::TYPE::data,
            "[UniformSegmentSetSampler] Undefined Extended Bounding "
            "Box (volume ==0)." );
        this->set_log_pdf( -std::log( volume ) );
    }

    OwnerSegment2D UniformSegmentSetSampler::sample(
        RandomEngine& engine ) const
    {
        auto seg = SegmentUniformSampler::sample(
            engine, domain_.extended_box(), length_, azimuth_ );
        return seg;
    }

    // Each change move leaves the reference measure invariant (uniform anchor
    // x length distribution x azimuth distribution), so the MH acceptance
    // only depends on the energy variation. The move type is drawn
    // independently of the current object.
    OwnerSegment2D UniformSegmentSetSampler::change(
        const OwnerSegment2D& obj, RandomEngine& engine ) const
    {
        if( engine.sample_bernoulli( translation_probability_ ) )
        {
            return translate( obj, engine );
        }
        if( engine.sample_bernoulli( rotation_probability_ ) )
        {
            return rotate( obj, engine );
        }
        return stretch( obj, engine );
    }

    OwnerSegment2D UniformSegmentSetSampler::translate(
        const OwnerSegment2D& obj, RandomEngine& engine ) const
    {
        // Symmetric move: the ball radius only depends on the length, which is
        // invariant by translation.
        const auto& extremities = obj.vertices();
        const geode::Sphere< 2 > ball{ extremities[0],
            move_ratio_ * obj.length() };
        const auto new_anchor =
            PointUniformSampler::sample< 2 >( engine, ball );
        if( !domain_.extended_contains( new_anchor ) )
        {
            // Out of the sampling domain: equivalent to a rejected move.
            // Re-sampling instead would break the move symmetry.
            return obj;
        }
        const auto translation = new_anchor - extremities[0];
        return OwnerSegment2D{ new_anchor, extremities[1] + translation };
    }

    OwnerSegment2D UniformSegmentSetSampler::rotate(
        const OwnerSegment2D& obj, RandomEngine& engine ) const
    {
        // Independent proposal of the azimuth from its distribution: the
        // proposal density cancels with the reference measure density.
        const auto& anchor = obj.vertices()[0];
        const auto azimuth = DoubleSampler::sample( engine, azimuth_ );
        return OwnerSegment2D{ anchor,
            anchor + direction_from_azimuth( azimuth ) * obj.length() };
    }

    OwnerSegment2D UniformSegmentSetSampler::stretch(
        const OwnerSegment2D& obj, RandomEngine& engine ) const
    {
        // Independent proposal of the length from its distribution: the
        // proposal density cancels with the reference measure density.
        const auto& anchor = obj.vertices()[0];
        const auto length = DoubleSampler::sample( engine, length_ );
        return OwnerSegment2D{ anchor,
            anchor + obj.normalized_direction() * length };
    }

    bool UniformSegmentSetSampler::is_valid_object(
        const OwnerSegment2D& obj ) const
    {
        // a segment belongs to the process if its anchor is in the extended
        // domain (same convention as the sampling and the change moves)
        return domain_.extended_contains( obj.vertices()[0] );
    }

    template <>
    std::unique_ptr< ObjectSetSampler< OwnerSegment2D > >
        build_objectset_sampler( const SpatialDomain< 2 >& domain,
            const ObjectSamplerConfig< OwnerSegment2D >& config )
    {
        return std::make_unique< UniformSegmentSetSampler >( domain, config );
    }

} // namespace geode