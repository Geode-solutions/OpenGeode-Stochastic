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
#include <numeric>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/direct/bounding_box_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/distance.hpp>

const int NUMBER_OF_SAMPLES = 10000;

template < geode::index_t dimension >
void test_sample_bounding_box(
    geode::RandomEngine& engine, const geode::BoundingBox< dimension >& box )
{
    geode::BoundingBoxSampler< dimension > spec_box{ box };

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value = spec_box.sample_uniform( engine );
        geode::OpenGeodeStochasticStochasticException::test(
            box.contains( value ), "[Point Box sampler] - point out of box." );
    }
}

void test_box_sampling1D()
{
    geode::RandomEngine random_engine;

    geode::Point1D min_point{ { 0. } };
    geode::Point1D max_point{ { 10. } };

    geode::BoundingBox1D box;
    box.add_point( min_point );
    box.add_point( max_point );
    test_sample_bounding_box< 1 >( random_engine, box );
}

void test_box_sampling2D()
{
    geode::RandomEngine random_engine;

    geode::Point2D min_point{ { 0., 0. } };
    geode::Point2D max_point{ { 10., 100. } };

    geode::BoundingBox2D box;
    box.add_point( min_point );
    box.add_point( max_point );
    test_sample_bounding_box< 2 >( random_engine, box );
}

void test_box_sampling3D()
{
    geode::RandomEngine random_engine;

    geode::Point3D min_point{ { 0., 0., 0. } };
    geode::Point3D max_point{ { 10., 100., -100. } };

    geode::BoundingBox3D box;
    box.add_point( min_point );
    box.add_point( max_point );
    test_sample_bounding_box< 3 >( random_engine, box );
}

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();

        test_box_sampling1D();
        test_box_sampling2D();
        test_box_sampling3D();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
