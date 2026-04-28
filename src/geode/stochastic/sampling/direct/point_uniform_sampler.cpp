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

#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/direct/ball_sampler.hpp>
#include <geode/stochastic/sampling/direct/bounding_box_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{

    template < index_t dimension >
    Point< dimension > PointUniformSampler::sample(
        RandomEngine& engine, const Object< dimension >& object )
    {
        return std::visit(
            [&engine]( auto&& obj ) {
                using OBJ = std::decay_t< decltype( obj ) >;
                if constexpr( std::is_same_v< OBJ, BoundingBox< dimension > > )
                {
                    BoundingBoxSampler< dimension > sampler{ obj };
                    return sampler.sample_uniform( engine );
                }
                if constexpr( std::is_same_v< OBJ, Sphere< dimension > > )
                {
                    BallSampler< dimension > sampler{ obj };
                    return sampler.sample_uniform( engine );
                }
                throw OpenGeodeStochasticStochasticException{ nullptr,
                    OpenGeodeException::TYPE::data,
                    "PointUniformSampler - Unsupported object for point "
                    "sampling" };
            },
            object );
    }

    template Point< 2 >
        opengeode_stochastic_stochastic_api PointUniformSampler::sample< 2 >(
            RandomEngine& engine, const Object< 2 >& object );
    template Point< 3 >
        opengeode_stochastic_stochastic_api PointUniformSampler::sample< 3 >(
            RandomEngine& engine, const Object< 3 >& object );

} // namespace geode