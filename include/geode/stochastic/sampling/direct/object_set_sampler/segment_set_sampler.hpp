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

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/basic_objects/sphere.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>
#include <geode/stochastic/sampling/direct/segment_uniform_sampler.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    class UniformSegmentSetSampler : public ObjectSetSampler< OwnerSegment2D >
    {
    public:
        UniformSegmentSetSampler( const SpatialDomain< 2 >& domain,
            const DoubleSampler::Distribution& length,
            const DoubleSampler::Distribution& azimuth )
            : ObjectSetSampler< OwnerSegment2D >{},
              domain_{ domain },
              length_{ length },
              azimuth_{ azimuth }
        {
            auto volume = domain_.extended_n_volume();
            OpenGeodeStochasticStochasticException::check_exception(
                volume != 0., nullptr, OpenGeodeException::TYPE::data,
                "[SegmentSetSampler] - Undefined Extended Bounding "
                "Box (volume ==0)." );
            this->log_pdf_ = -std::log( volume );
        }

        OwnerSegment2D sample( RandomEngine& engine ) const override
        {
            auto seg = SegmentUniformSampler::sample(
                engine, domain_.extended_box(), length_, azimuth_ );
            return seg;
        }

        OwnerSegment2D change(
            const OwnerSegment2D& obj, RandomEngine& engine ) const override
        {
            double ratio = 0.1;
            const auto& extremities = obj.vertices();
            const auto current =
                static_cast< local_index_t >( engine.sample_bernoulli( 0.5 ) );
            const auto other = 1 - current;

            geode::Sphere< 2 > ball{ extremities[current],
                ratio * obj.length() };

            auto new_point = PointUniformSampler::sample< 2 >( engine, ball );
            constexpr index_t max_try{ 100 };
            for( const auto try_id : geode::Range{ max_try } )
            {
                if( domain_.extended_contains( new_point )
                    || domain_.extended_contains( extremities[other] ) )
                {
                    OwnerSegment2D new_segment{ obj };
                    new_segment.set_point( current, new_point );
                    return new_segment;
                }
                new_point = PointUniformSampler::sample< 2 >( engine, ball );
            }
            throw OpenGeodeStochasticStochasticException{ nullptr,
                OpenGeodeException::TYPE::internal,
                "[SegmentSetSampler] - Cannot find a point in the box" };
            return obj;
        }

    private:
        bool is_valid_object( const OwnerSegment2D& obj ) const override
        {
            const auto& extremities = obj.vertices();
            return domain_.extended_contains( extremities[0] )
                   || domain_.extended_contains( extremities[1] );
        }

    private:
        const SpatialDomain< 2 >& domain_;
        DoubleSampler::Distribution length_;
        DoubleSampler::Distribution azimuth_;
    };

} // namespace geode