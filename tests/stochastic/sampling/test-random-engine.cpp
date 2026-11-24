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
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <optional>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/distributions.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>
// NOLINTBEGIN(*-magic-numbers)
namespace
{
    const int NUMBER_OF_SAMPLES = 100000;

    void test_reproducibility()
    {
        std::vector< std::string > seeds = { "same-seed", "geode-solutions",
            "dfgzejhèçsodj", "&%;:!§" };
        for( const auto& seed : seeds )
        { // Create first engine geode::
            geode::RandomEngine engine1;
            engine1.set_seed( seed );

            // Create second engine with same seed
            geode::RandomEngine engine2;
            engine2.set_seed( seed );

            // Define a uniform distribution
            geode::UniformClosed< int > dist;
            dist.min_value = 1;
            dist.max_value = 100;

            // Sample more values to check sequence reproducibility
            for( const auto value : geode::Range{ 100 } )
            {
                geode_unused( value );
                geode::OpenGeodeStochasticStochasticException::test(
                    engine1.sample_uniform( dist )
                        == engine2.sample_uniform( dist ),
                    "[REPRODUCIBILITY] - Same seed should produce same "
                    "output." );
            }
            geode::Logger::info( "Test Reproducibility: ", seed, " SUCCESS " );
        }
        geode::Logger::info( "Test Reproducibility: SUCCESS " );
    }

    // A textual seed is hashed with FNV-1a 64 bits, which is stable across
    // runs and platforms: seeding with the word or with its hash value must
    // give the same sequence. Reference values from the FNV-1a specification
    // ("", "a") and computed independently (others).
    void test_textual_seed_hash()
    {
        const std::vector< std::pair< std::string_view, uint64_t > >
            references = { { "", 0xcbf29ce484222325 },
                { "a", 0xaf63dc4c8601ec8c },
                { "same-seed", 0x7b87eb74217c1663 },
                { "geode-solutions", 0x3ce1c85d8114d36c } };
        geode::UniformClosed< int > dist;
        dist.min_value = 0;
        dist.max_value = 1000000000;
        for( const auto& [word, hash] : references )
        {
            geode::RandomEngine word_engine;
            word_engine.set_seed( word );
            geode::RandomEngine hash_engine;
            hash_engine.set_seed( hash );
            for( const auto value : geode::Range{ 10 } )
            {
                geode_unused( value );
                geode::OpenGeodeStochasticStochasticException::test(
                    word_engine.sample_uniform( dist )
                        == hash_engine.sample_uniform( dist ),
                    "[TextualSeed] - Seed \"", word,
                    "\" should be hashed to its FNV-1a 64 bits value." );
            }
        }
        geode::Logger::info( "Test textual seed hash: SUCCESS " );
    }

    // A seeded engine must give the same sequence from one run to another:
    // the expected values are hard-coded. They come from std::mt19937_64
    // seeded by std::seed_seq (both specified by the C++ standard) and
    // absl::Uniform, computed outside of RandomEngine.
    void test_cross_run_reproducibility()
    {
        geode::UniformClosed< int > dist;
        dist.min_value = 0;
        dist.max_value = 1000000000;
        const auto check = []( geode::RandomEngine& engine,
                               const geode::UniformClosed< int >& law,
                               const std::array< int, 5 >& expected,
                               std::string_view seed ) {
            for( const auto value : expected )
            {
                const auto sample = engine.sample_uniform( law );
                geode::OpenGeodeStochasticStochasticException::test(
                    sample == value, "[CrossRunReproducibility] - Seed ", seed,
                    ": expected ", value, ", got ", sample );
            }
        };
        geode::RandomEngine number_engine;
        number_engine.set_seed( uint64_t{ 42 } );
        check( number_engine, dist,
            { 213580860, 705506416, 727284205, 46917137, 849020427 }, "42" );
        geode::RandomEngine word_engine;
        word_engine.set_seed( "same-seed" );
        check( word_engine, dist,
            { 487694039, 815218771, 864671017, 660062510, 339159719 },
            "same-seed" );
        geode::Logger::info( "Test cross run reproducibility: SUCCESS " );
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
        for( auto& x : data )
        {
            double diff = x - mean;
            accum += diff * diff;
        }
        return accum / ( data.size() - 1 );
    }
    // NOLINTBEGING(*-magic-numbers)
    template < typename T >
    void test_distribution_mean_and_variance( const std::vector< T >& data,
        double expected_mean,
        double expected_var,
        double k_coef = 5.0 )
    {
        const auto data_size = data.size();
        const double mean = compute_mean( data );
        const double variance = compute_variance( data, mean );
        geode::Logger::info( "Test Distribution ",
            ": mean / expected mean = ", mean, "/", expected_mean,
            " variance / expected variance = ", variance, "/", expected_var );

        const double se_mean = std::sqrt( expected_var / data_size );
        const double se_var =
            std::sqrt( 2.0 * expected_var * expected_var / ( data_size - 1 ) );

        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( mean - expected_mean ) < k_coef * se_mean,
            "[Uniform] - Wrong expected mean." );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( variance - expected_var ) < k_coef * se_var,
            "[Uniform] - Wrong expected std." );
    }

    template < typename T >
    double compute_expected_variance( T min_value, T max_value )
    {
        if constexpr( std::is_integral_v< T > )
        {
            const auto distance =
                static_cast< double >( max_value - min_value + 1 );
            return ( distance * distance - 1.0 ) / 12.0;
        }
        else if constexpr( std::is_floating_point_v< T > )
        {
            const auto distance =
                static_cast< double >( max_value - min_value );
            return ( distance * distance ) / 12.0;
        }
        else
        {
            static_assert( std::is_arithmetic_v< T >,
                "Unsupported type to compute theoretical variance." );
        }
    }

    template < typename Type >
    void test_uniform( const Type min_value,
        const Type max_value,
        geode::RandomEngine& engine )
    {
        std::vector< Type > samples;
        samples.reserve( NUMBER_OF_SAMPLES );

        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            geode::UniformClosed< Type > dist_closed;
            dist_closed.min_value = min_value;
            dist_closed.max_value = max_value;

            Type value = engine.sample_uniform( dist_closed );
            samples.emplace_back( value );
            geode::OpenGeodeStochasticStochasticException::test(
                value >= min_value && value <= max_value,
                "[Uniform] -  value out of range." );
        }

        double expected_mean = ( min_value + max_value ) / 2.0;
        double expected_var = compute_expected_variance( min_value, max_value );

        test_distribution_mean_and_variance(
            samples, expected_mean, expected_var );
    }

    void test_gaussian(
        double mean_value, double std_value, geode::RandomEngine& engine )
    {
        std::vector< double > samples;
        samples.reserve( NUMBER_OF_SAMPLES );
        geode::Gaussian spec;
        spec.mean = mean_value;
        spec.standard_deviation = std_value;

        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            double value = engine.sample_gaussian( spec );
            samples.emplace_back( value );
        }

        test_distribution_mean_and_variance( samples, mean_value, std_value );
    }

    void test_truncated_gaussian( double mean_value,
        double std_value,
        std::optional< double > min_value,
        std::optional< double > max_value,
        geode::RandomEngine& engine )
    {
        const auto max =
            max_value.value_or( std::numeric_limits< double >::infinity() );
        const auto min =
            min_value.value_or( -std::numeric_limits< double >::infinity() );

        std::vector< double > samples;
        samples.reserve( NUMBER_OF_SAMPLES );

        geode::TruncatedGaussian spec;
        spec.mean = mean_value;
        spec.standard_deviation = std_value;
        spec.min_value = min_value;
        spec.max_value = max_value;

        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            double value = engine.sample_truncated_gaussian( spec );
            samples.emplace_back( value );
            geode::OpenGeodeStochasticStochasticException::test(
                value >= min && value <= max,
                "[Gaussian] -  value out of range." );
        }

        // Implementer un test de Chi2 ou Kolmogorov–Smirnov
    }

    void test_bernoulli(
        double probability_of_success, geode::RandomEngine& engine )
    {
        int success_count = 0;

        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            if( engine.sample_bernoulli( probability_of_success ) )
            {
                ++success_count;
            }
        }

        const double empirical_probability =
            static_cast< double >( success_count ) / NUMBER_OF_SAMPLES;

        geode::Logger::info( "Test Bernoulli ",
            ": empirical / expected = ", empirical_probability, " / ",
            probability_of_success );

        // NOLINTBEGING(*-magic-numbers)
        geode::OpenGeodeStochasticStochasticException::test(
            abs( empirical_probability - probability_of_success ) < 0.05,
            "[Bernoulli] - Empirical probability out of tolerance." );
    }

    void test_gaussian_tail( double lower, geode::RandomEngine& engine )
    {
        // N(0,1) truncated to [lower, inf): mean = phi(lower) / Q(lower)
        geode::TruncatedGaussian law;
        law.min_value = lower;
        double sum{ 0. };
        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            sum += engine.sample_truncated_gaussian( law );
        }
        const auto density =
            std::exp( -0.5 * lower * lower ) / std::sqrt( 2. * M_PI );
        const auto tail = 0.5 * std::erfc( lower / std::sqrt( 2. ) );
        const auto expected = density / tail;
        const auto mean = sum / NUMBER_OF_SAMPLES;
        geode::Logger::info( "Test Gaussian tail above ", lower,
            ": empirical / expected mean = ", mean, " / ", expected );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( mean - expected ) < 0.01,
            "[Gaussian] - Tail mean out of tolerance." );
    }

    template < typename Sample >
    bool throws_data_error( const Sample& sample )
    {
        try
        {
            sample();
        }
        catch( const geode::OpenGeodeException& exception )
        {
            return exception.type() == geode::OpenGeodeException::TYPE::data;
        }
        return false;
    }

    template < typename Sample >
    bool throws_within_draws( const Sample& sample )
    {
        for( const auto draw : geode::Range{ 100 } )
        {
            geode_unused( draw );
            if( throws_data_error( sample ) )
            {
                return true;
            }
        }
        return false;
    }

    void test_unrepresentable_samples_are_rejected(
        geode::RandomEngine& engine )
    {
        // Valid laws whose values can exceed the double range: an error is
        // raised instead of returning inf or NaN.
        geode::TruncatedPowerLaw power_law;
        power_law.alpha = 1. + 1e-7;
        power_law.min_value = 1.;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_within_draws( [&] {
                return engine.sample_truncated_powerlaw( power_law );
            } ),
            "[PowerLaw] - Infinite sample returned." );
        geode::Gaussian gaussian;
        gaussian.mean = 1e308;
        gaussian.standard_deviation = 1e308;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_within_draws( [&] {
                return engine.sample_gaussian( gaussian );
            } ),
            "[Gaussian] - Infinite sample returned." );
        geode::TruncatedGaussian truncated;
        truncated.mean = 1e308;
        truncated.standard_deviation = 1e308;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_within_draws( [&] {
                return engine.sample_truncated_gaussian( truncated );
            } ),
            "[TruncatedGaussian] - Infinite sample returned." );
        geode::TruncatedLogNormal lognormal;
        lognormal.mean = 700.;
        lognormal.standard_deviation = 10.;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_within_draws( [&] {
                return engine.sample_truncated_lognormal( lognormal );
            } ),
            "[LogNormal] - Infinite sample returned." );
    }

    void test_lognormal_tail( geode::RandomEngine& engine )
    {
        // log( x ) truncated 20 standard deviations above its mean: the
        // samples are correct; far beyond, an explicit error is raised.
        geode::TruncatedLogNormal law;
        law.min_value = std::exp( 20. );
        for( const auto count : geode::Range{ 1000 } )
        {
            geode_unused( count );
            const auto value = engine.sample_truncated_lognormal( law );
            geode::OpenGeodeStochasticStochasticException::test(
                value >= std::exp( 20. ) && value < std::exp( 21. ),
                "[LogNormal] - Tail sample out of range: ", value );
        }
        law.min_value = 1e20;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_data_error( [&] {
                return engine.sample_truncated_lognormal( law );
            } ),
            "[LogNormal] - No error for a range 46 sigma above the mean." );
    }

    void test_invalid_laws_are_rejected( geode::RandomEngine& engine )
    {
        // The engine applies the rules of is_valid(): an invalid law is
        // rejected instead of giving out of range samples.
        geode::TruncatedPowerLaw power_law;
        power_law.min_value = -1.;
        power_law.max_value = 2.;
        geode::OpenGeodeStochasticStochasticException::test(
            throws_data_error( [&] {
                return engine.sample_truncated_powerlaw( power_law );
            } ),
            "[PowerLaw] - Negative lower bound accepted." );
        power_law.min_value = 1.;
        power_law.alpha = std::numeric_limits< double >::infinity();
        geode::OpenGeodeStochasticStochasticException::test(
            throws_data_error( [&] {
                return engine.sample_truncated_powerlaw( power_law );
            } ),
            "[PowerLaw] - Infinite exponent accepted." );
        // x^-alpha with alpha <= 1 cannot be normalized up to infinity
        for( const auto alpha : { 1., 0.5 } )
        {
            power_law.alpha = alpha;
            power_law.max_value = std::nullopt;
            geode::OpenGeodeStochasticStochasticException::test(
                throws_data_error( [&] {
                    return engine.sample_truncated_powerlaw( power_law );
                } ),
                "[PowerLaw] - Exponent ", alpha,
                " accepted without upper bound." );
            geode::OpenGeodeStochasticStochasticException::test(
                !power_law.is_valid(), "[PowerLaw] - Exponent ", alpha,
                " is valid without upper bound." );
        }
    }

    void test_truncated_powerlaw( double alpha,
        double min_value,
        std::optional< double > max_value,
        geode::RandomEngine& engine )
    {
        geode::TruncatedPowerLaw law;
        law.alpha = alpha;
        law.min_value = min_value;
        law.max_value = max_value;
        const auto max =
            max_value.value_or( std::numeric_limits< double >::infinity() );
        std::vector< double > samples;
        samples.reserve( NUMBER_OF_SAMPLES );
        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            const auto value = engine.sample_truncated_powerlaw( law );
            geode::OpenGeodeStochasticStochasticException::test(
                value >= min_value && value <= max, "[PowerLaw] - value ",
                value, " out of [", min_value, ",", max, "]." );
            samples.push_back( value );
        }
        // median of the density proportional to x^-alpha on [min, max]
        const auto median = [&] {
            if( alpha == 1. )
            {
                return std::sqrt( min_value * max );
            }
            // with the ratio max / min: min^(1-alpha) may overflow
            const auto ratio = std::pow( max / min_value, 1. - alpha );
            return min_value
                   * std::pow( 0.5 * ( 1. + ratio ), 1. / ( 1. - alpha ) );
        }();
        std::nth_element( samples.begin(),
            samples.begin() + NUMBER_OF_SAMPLES / 2, samples.end() );
        const auto empirical_median = samples[NUMBER_OF_SAMPLES / 2];
        geode::Logger::info( "Test PowerLaw alpha=", alpha,
            ": empirical / expected median = ", empirical_median, " / ",
            median );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( empirical_median - median ) < 0.02 * median,
            "[PowerLaw] - Empirical median out of tolerance." );
    }

    // Samples are wrapped in [0, 2pi[ whatever the mean, also a negative one
    // or one beyond 2pi, and are centered on the mean direction.
    void test_von_mises(
        double mean, double kappa, geode::RandomEngine& engine )
    {
        geode::VonMises law;
        law.mean = mean;
        law.concentration = kappa;
        constexpr geode::index_t NUMBER_OF_SAMPLES{ 5000 };
        double sum_sin{ 0. };
        double sum_cos{ 0. };
        for( const auto count : geode::Range{ NUMBER_OF_SAMPLES } )
        {
            geode_unused( count );
            const auto value = engine.sample_von_mises( law );
            geode::OpenGeodeStochasticStochasticException::test(
                value >= 0. && value < 2. * M_PI, "[VonMises] - mean ", mean,
                ", kappa ", kappa, ": sample ", value, " out of [0, 2pi[." );
            sum_sin += std::sin( value );
            sum_cos += std::cos( value );
        }
        const auto difference = std::atan2( sum_sin, sum_cos ) - mean;
        const auto gap = std::abs(
            std::atan2( std::sin( difference ), std::cos( difference ) ) );
        geode::OpenGeodeStochasticStochasticException::test( gap < 0.1,
            "[VonMises] - mean ", mean, ", kappa ", kappa,
            ": samples are centered ", gap, " away from the mean direction." );
    }
} // namespace
int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::RandomEngine random_engine;
        test_reproducibility();
        test_textual_seed_hash();
        test_cross_run_reproducibility();

        random_engine.set_seed( "@-test-radom-engine@" );

        geode::Logger::info( "TEST UNIFORM SAMPLING" );
        test_uniform< geode::index_t >( 1, 15, random_engine );
        test_uniform< geode::local_index_t >( 1, 6, random_engine );
        test_uniform< geode::signed_index_t >( -18, 6, random_engine );
        test_uniform< float >( 15.545, 18.9524, random_engine );
        test_uniform< double >( -100.54, 100.6, random_engine );
        geode::Logger::info( "TEST UNIFORM SAMPLING - SUCCESS" );

        geode::Logger::info( "TEST GAUSSIAN SAMPLING" );
        test_gaussian( 0., 1., random_engine );
        test_gaussian( -100., 1., random_engine );
        geode::Logger::info( "TEST GAUSSIAN SAMPLING - SUCCESS" );

        geode::Logger::info( "TEST TRUNCATED GAUSSIAN SAMPLING" );
        test_truncated_gaussian(
            0., 1., std::nullopt, std::nullopt, random_engine );
        test_truncated_gaussian( 0., 1., std::nullopt, 1., random_engine );
        test_truncated_gaussian( 0., 1., -1., std::nullopt, random_engine );
        test_truncated_gaussian( 10., 1., 5., 15., random_engine );
        // tails up to 6 standard deviations are valid laws (is_valid)
        test_truncated_gaussian( 0., 1., 5.5, std::nullopt, random_engine );
        test_truncated_gaussian( 0., 1., std::nullopt, -5.5, random_engine );
        test_truncated_gaussian( 0., 1., 6., 8., random_engine );
        test_gaussian_tail( 5.5, random_engine );
        geode::Logger::info( "TEST TRUNCATED GAUSSIAN SAMPLING - SUCCESS" );

        geode::Logger::info( "TEST POWER LAW SAMPLING" );
        test_truncated_powerlaw( 2., 1., 10., random_engine );
        test_truncated_powerlaw( 2., 1., std::nullopt, random_engine );
        test_truncated_powerlaw( 3.5, 0.5, std::nullopt, random_engine );
        test_truncated_powerlaw( 1., 1., 10., random_engine );
        test_truncated_powerlaw( 0.5, 1., 10., random_engine );
        // a large exponent with a small lower bound: xmin^(1-alpha) overflows
        test_truncated_powerlaw( 500., 1e-6, std::nullopt, random_engine );
        test_truncated_powerlaw( 500., 1e-6, 1., random_engine );
        test_invalid_laws_are_rejected( random_engine );
        test_unrepresentable_samples_are_rejected( random_engine );
        test_lognormal_tail( random_engine );
        geode::Logger::info( "TEST POWER LAW SAMPLING - SUCCESS" );

        geode::Logger::info( "TEST VON MISES SAMPLING" );
        for( const auto mean : { -10., -4., 0., 1., 7., 100. } )
        {
            test_von_mises( mean, 1., random_engine );
            test_von_mises( mean, 5000., random_engine );
        }
        geode::Logger::info( "TEST VON MISES SAMPLING - SUCCESS" );

        geode::Logger::info( "TEST BERNOUILLI SAMPLING" );
        test_bernoulli( 0.5, random_engine );
        test_bernoulli( 0.1, random_engine );
        test_bernoulli( 0.25, random_engine );
        test_bernoulli( 0.06, random_engine );
        test_bernoulli( 0.96, random_engine );
        geode::Logger::info( "TEST BERNOUILLI SAMPLING - SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
// NOLINTEND(*-magic-numbers)
