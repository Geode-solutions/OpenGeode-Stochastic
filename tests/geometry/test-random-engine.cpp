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

#include <geode/stochastic/geometry/random_engine.hpp>
#include <geode/stochastic/geometry/random_sample_spec.hpp>

template < typename T >
double compute_mean( const std::vector< T >& data )
{
    double sum = std::accumulate( data.begin(), data.end(), 0.0 );
    return sum / data.size();
}

template < typename T >
double compute_variance( const std::vector< T >& data, double mean )
{
    double accum = 0.0;
    for( auto x : data )
    {
        double diff = x - mean;
        accum += diff * diff;
    }
    return accum / ( data.size() - 1 );
}

void test_uniform_int( geode::RandomEngine& engine )
{
    const int number_of_samples = 10000;
    const int min_value = 1;
    const int max_value = 6;
    std::vector< int > samples{ number_of_samples };
    geode::UniformInt spec{ min_value, max_value };

    for( const auto i : geode::Range{ number_of_samples } )
    {
        int value = engine.sample( geode::UniformInt{ min_value, max_value } );
        samples.emplace_back( value );
        OPENGEODE_ASSERT( value >= min_value && value <= max_value,
            "[UniformInt] -  value out of range." );
    }
    double mean = compute_mean( samples );
    double variance = compute_variance( samples, mean );

    double expected_mean = ( min_value + max_value ) / 2.0;
    double expected_var =
        ( ( max_value - min_value + 1 ) * ( max_value - min_value + 1 ) - 1 )
        / 12.0;
    OPENGEODE_ASSERT( mean > expected_mean - 0.1 && mean < expected_mean + 0.1,
        "[UniformInt] - Wrong expected mean." );
    OPENGEODE_ASSERT(
        variance > expected_var - 0.1 && variance < expected_var + 0.1,
        "[UniformInt] - Wrong expected std." );
    geode::Logger::info( "Test UniformInt: SUCCESS" );
}

int main()
{
    try
    {
        geode::StochasticGeometryLibrary::initialize();
        geode::RandomEngine random_engine( 27041986 );
        test_uniform_int( random_engine );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
