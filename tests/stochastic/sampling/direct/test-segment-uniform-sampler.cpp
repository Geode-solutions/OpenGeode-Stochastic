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
//#include <numeric>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>
#include <geode/stochastic/sampling/direct/segment_uniform_sampler.hpp>
#include <geode/stochastic/sampling/distributions.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>

const int NUMBER_OF_SAMPLES = 10000;

void test_sample_segment(
    geode::RandomEngine& engine, const geode::BoundingBox< 2 >& box )
{
    geode::UniformClosed< double > length;
    geode::UniformClosed< double > az;

    for( const auto i : geode::Range{ NUMBER_OF_SAMPLES } )
    {
        auto value =
            geode::SegmentUniformSampler::sample( engine, box, length, az );
        // OPENGEODE_ASSERT( box.contains( value ),
        //     "[SegmentUniformSampler] - segment out of box." );
    }
}
void test_segment_sampling2D()
{
    geode::RandomEngine random_engine;

    geode::Point2D min_point{ { 0.1, 0.8 } };
    geode::Point2D max_point{ { 1.2, 1. } };
    geode::BoundingBox2D box;
    box.add_point( min_point );
    box.add_point( max_point );
    test_sample_segment( random_engine, box );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        test_segment_sampling2D();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
