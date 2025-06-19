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

void test_reproducibility()
{
    std::vector< std::string > seeds = { "same-seed", "geode-solutions",
        "dfgzejhèçsodj", "&%;:!§" };
    for( const auto seed : seeds )
    { // Create first engine geode::
        geode::RandomEngine engine1;
        engine1.set_seed( seed );

        // Create second engine with same seed
        geode::RandomEngine engine2;
        engine2.set_seed( seed );

        // Define a uniform distribution
        geode::Uniform< int > dist{ 1, 100 };

        // Sample more values to check sequence reproducibility
        for( const auto value : geode::Range{ 100 } )
        {
            OPENGEODE_ASSERT( engine1.sample_uniform( dist )
                                  == engine2.sample_uniform( dist ),
                "[REPRODUCIBILITY] - Same seed should produce same output." );
        }
        geode::Logger::info( "Test Reproducibility: ", seed, " SUCCESS " );
    }
    geode::Logger::info( "Test Reproducibility: SUCCESS " );
}

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
template < typename Type >
void test_uniform(
    const Type min_value, const Type max_value, geode::RandomEngine& engine )
{
    const int number_of_samples = 10000;
    std::vector< Type > samples;
    samples.reserve( number_of_samples );

    for( const auto i : geode::Range{ number_of_samples } )
    {
        Type value = engine.sample_uniform(
            geode::Uniform< Type >{ min_value, max_value } );
        samples.emplace_back( value );
        OPENGEODE_ASSERT( value >= min_value && value <= max_value,
            "[Uniform] -  value out of range." );
    }
    double mean = compute_mean( samples );
    double variance = compute_variance( samples, mean );

    double expected_mean = ( min_value + max_value ) / 2.0;
    double expected_var =
        ( ( max_value - min_value + 1 ) * ( max_value - min_value + 1 ) - 1 )
        / 12.0;
    OPENGEODE_ASSERT( mean > expected_mean - 0.1 && mean < expected_mean + 0.1,
        "[Uniform] - Wrong expected mean." );
    OPENGEODE_ASSERT(
        variance > expected_var - 0.1 && variance < expected_var + 0.1,
        "[Uniform] - Wrong expected std." );
    geode::Logger::info( "Test Uniform ",
        ": SUCCESS - mean / expected mean = ", mean, "/", expected_mean,
        " variance / expected variance = ", variance, "/", expected_var );
}

void test_gaussian( double mean_value,
    double std_value,
    std::optional< double > min_value,
    std::optional< double > max_value,
    geode::RandomEngine& engine )
{
    const auto max =
        max_value.value_or( std::numeric_limits< double >::infinity() );
    const auto min =
        min_value.value_or( -std::numeric_limits< double >::infinity() );
    OPENGEODE_ASSERT( 2 * mean_value - min - max < geode::GLOBAL_EPSILON,
        "[Gaussian] -  Test not relevant min an max are not symetric "
        "compare to mean value." );

    const int number_of_samples = 10000;
    std::vector< double > samples;
    samples.reserve( number_of_samples );
    geode::Gaussian spec{ mean_value, std_value };
    spec.min = min_value;
    spec.max = max_value;

    for( const auto i : geode::Range{ number_of_samples } )
    {
        double value = engine.sample_gaussian( spec );
        samples.emplace_back( value );
        OPENGEODE_ASSERT(
            value >= min && value <= max, "[Gaussian] -  value out of range." );
    }
    double mean = compute_mean( samples );
    double variance = compute_variance( samples, mean );

    double expected_mean = mean_value;
    double expected_var = std_value * std_value;

    OPENGEODE_ASSERT( mean > expected_mean - 0.1 && mean < expected_mean + 0.1,
        "[Gaussian] - Wrong expected mean." );
    OPENGEODE_ASSERT(
        variance > expected_var - 0.1 && variance < expected_var + 0.1,
        "[Gaussian] - Wrong expected std." );
    geode::Logger::info( "Test Gaussian ",
        ": SUCCESS - mean / expected mean = ", mean, "/", expected_mean,
        " variance / expected variance = ", variance, "/", expected_var );
}

void test_bernoulli(
    double probability_of_success, geode::RandomEngine& engine )
{
    const int number_of_samples = 10000;
    int success_count = 0;

    for( const auto i : geode::Range{ number_of_samples } )
    {
        if( engine.sample_bernoulli( probability_of_success ) )
        {
            ++success_count;
        }
    }

    const double empirical_probability =
        static_cast< double >( success_count ) / number_of_samples;

    OPENGEODE_ASSERT(
        abs( empirical_probability - probability_of_success ) < 0.05,
        "[Bernoulli] - Empirical probability out of tolerance." );

    geode::Logger::info( "Test Bernoulli ",
        ": SUCCESS - empirical / expected = ", empirical_probability, " / ",
        probability_of_success );
}

int main()
{
    try
    {
        geode::StochasticGeometryLibrary::initialize();
        geode::RandomEngine random_engine;
        test_reproducibility();

        test_uniform< geode::index_t >( 1, 15, random_engine );
        test_uniform< geode::local_index_t >( 1, 6, random_engine );
        test_uniform< geode::signed_index_t >( -18, 6, random_engine );
        test_uniform< float >( 15.545, 18.9524, random_engine );
        test_uniform< double >( -100.54, 100.6, random_engine );

        test_gaussian( 0., 1., std::nullopt, std::nullopt, random_engine );
        test_gaussian( 10., 1., 5., 15., random_engine );

        test_bernoulli( 0.5, random_engine );
        test_bernoulli( 0.1, random_engine );
        test_bernoulli( 0.25, random_engine );
        test_bernoulli( 0.06, random_engine );
        test_bernoulli( 0.96, random_engine );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
