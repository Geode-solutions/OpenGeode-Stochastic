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

#include <geode/stochastic/sampling/direct/segment_uniform_sampler.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/angle.hpp>
#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/vector.hpp>

namespace
{
    geode::Vector2D direction_from_azimuth_angle( const geode::Angle &azimuth )
    {
        return geode::Vector2D( { azimuth.sin(), azimuth.cos() } );
    }
} // namespace
namespace geode
{

    OwnerSegment2D SegmentUniformSampler::sample( RandomEngine &engine,
        const PointUniformSampler::Object< 2 > &object,
        const DoubleSampler::Distribution &length,
        const DoubleSampler::Distribution &azimuth_rad )
    {
        auto point1 = PointUniformSampler::sample( engine, object );

        auto segment_length = DoubleSampler::sample( engine, length );
        auto segment_azimuth = Angle::create_from_radians(
            DoubleSampler::sample( engine, azimuth_rad ) );

        auto point2 =
            point1
            + direction_from_azimuth_angle( segment_azimuth ) * segment_length;

        return OwnerSegment2D{ point1, point2 };
    }

} // namespace geode