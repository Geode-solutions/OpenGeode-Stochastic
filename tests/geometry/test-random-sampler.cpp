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

#include <geode/stochastic/geometry/distributions.hpp>
#include <geode/stochastic/geometry/random_engine.hpp>
#include <geode/stochastic/geometry/random_sampler.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/distance.hpp>

const int NUMBER_OF_SAMPLES = 10000;

void test_double_sampler( geode::RandomEngine& engine )
{
    geode::UniformClosed< double > uniform_closed_double;
    double value = 1000.;
    value =
        geode::RandomSampler< double >::sample( engine, uniform_closed_double );
    OPENGEODE_ASSERT(
        value >= 0. && value <= 1., "[Uniform] -  value out of range." );

    geode::UniformClosedOpen< double > uniform_closedopen_double;
    value = 1000.;
    value = geode::RandomSampler< double >::sample(
        engine, uniform_closedopen_double );
    OPENGEODE_ASSERT(
        value >= 0. && value <= 1., "[Uniform] -  value out of range." );

    geode::TruncatedGaussian t_gaussian_double;
    t_gaussian_double.mean = 0.;
    t_gaussian_double.mean = 1.;
    t_gaussian_double.max_value = 1.;
    t_gaussian_double.min_value = -1.;
    value = 1000.;
    value = geode::RandomSampler< double >::sample( engine, t_gaussian_double );
    OPENGEODE_ASSERT( value >= -1. && value <= 1.,
        "[TruncatedGaussian] -  value out of range." );

    geode::Gaussian gaussian_double;
    t_gaussian_double.mean = 0.;
    t_gaussian_double.mean = 1.;

    value = 1000.;
    value = geode::RandomSampler< double >::sample( engine, t_gaussian_double );
    OPENGEODE_ASSERT( value != 1000., "[Gaussian] -  value did not changed." );
}

void test_point2D_sampler( geode::RandomEngine& engine )
{
    geode::Point2D min_point{ { 0., 0. } };
    geode::Point2D max_point{ { 100., -100. } };

    geode::BoundingBox2D box;
    box.add_point( min_point );
    box.add_point( max_point );

    geode::UniformBox< 2 > spec_box{ box };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value =
            geode::RandomSampler< geode::Point2D >::sample( engine, spec_box );
        OPENGEODE_ASSERT( box.contains( value ),
            "[Point Box sampler 2D] - point out of box." );
    }
    geode::UniformBall< 2 > spec_ball{ min_point, 1 };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value =
            geode::RandomSampler< geode::Point2D >::sample( engine, spec_ball );
        OPENGEODE_ASSERT( geode::point_point_distance( value, min_point ),
            "[Point Ball sampler 2D] - point too far from center." );
    }
}

void test_point3D_sampler( geode::RandomEngine& engine )
{
    geode::Point3D min_point{ { 0., 0., 0. } };
    geode::Point3D max_point{ { 10., 100., -100. } };

    geode::BoundingBox3D box;
    box.add_point( min_point );
    box.add_point( max_point );

    geode::UniformBox< 3 > spec_box{ box };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value =
            geode::RandomSampler< geode::Point3D >::sample( engine, spec_box );
        OPENGEODE_ASSERT( box.contains( value ),
            "[Point Box sampler 3D] - point out of box." );
    }
    geode::UniformBall< 3 > spec_ball{ min_point, 1 };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value =
            geode::RandomSampler< geode::Point3D >::sample( engine, spec_ball );
        OPENGEODE_ASSERT( geode::point_point_distance( value, min_point ),
            "[Point Ball sampler 3D] - point too far from center." );
    }
}

int main()
{
    try
    {
        geode::StochasticGeometryLibrary::initialize();
        geode::RandomEngine random_engine;

        test_double_sampler( random_engine );
        test_point2D_sampler( random_engine );
        test_point3D_sampler( random_engine );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
