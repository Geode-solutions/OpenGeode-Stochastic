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
#include <numeric>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/direct/ball_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/distance.hpp>

const int NUMBER_OF_SAMPLES = 10000;

template < geode::index_t dimension >
void test_sample_ball(
    geode::RandomEngine &engine, const geode::Sphere< dimension > &ball )
{
    geode::BallSampler< dimension > spec_ball{ ball };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value = spec_ball.sample_uniform( engine );
        OPENGEODE_EXCEPTION( geode::point_point_distance( value, ball.origin() )
                                 <= ball.radius() + geode::GLOBAL_EPSILON,
            "[Point Ball sampler] - point too far from center." );
    }
}

void test_ball_sampling2D()
{
    geode::RandomEngine random_engine;

    geode::Point2D center_point{ { 10., 100. } };
    geode::Sphere2D ball{ center_point, 0.5 };
    test_sample_ball< 2 >( random_engine, ball );
}

void test_ball_sampling3D()
{
    geode::RandomEngine random_engine;

    geode::Point3D center_point{ { 10., 100., -100. } };
    geode::Sphere3D ball{ center_point, 10 };
    test_sample_ball< 3 >( random_engine, ball );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        test_ball_sampling2D();
        test_ball_sampling3D();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
