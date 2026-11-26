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

#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/basic/pimpl_impl.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/types.hpp>
#include <geode/stochastic/common.hpp>

#include <absl/random/bernoulli_distribution.h>
#include <absl/random/gaussian_distribution.h>
#include <absl/random/random.h>
#include <absl/random/uniform_int_distribution.h>
#include <absl/random/uniform_real_distribution.h>

#include <algorithm>
#include <limits>
#include <random>

namespace
{
    // Wraps an angle in [0, 2pi[, whatever its sign or its number of turns.
    double wrap_to_two_pi( double angle )
    {
        constexpr auto TWO_PI = 2.0 * M_PI;
        auto wrapped = std::fmod( angle, TWO_PI );
        if( wrapped < 0.0 )
        {
            wrapped += TWO_PI;
            // a tiny negative remainder rounds to 2pi
            if( wrapped >= TWO_PI )
            {
                return 0.0;
            }
        }
        return wrapped;
    }

    double normal_cdf( double x )
    {
        return 0.5 * std::erfc( -x / std::sqrt( 2.0 ) );
    }

    // NormalQuantile based on Peter J. Acklam's inverse normal CDF
    // approximation. Original algorithm released to the public domain.
    // Reference:
    // https://web.archive.org/web/20150910005011/http://home.online.no/~pjacklam
    double normal_quantile( double p )
    {
        geode::OpenGeodeStochasticStochasticException::check_exception(
            p >= 0.0 && p <= 1.0, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[normal_quantile] p must be in (0,1). Check the consistencies "
            "between min,mean and max values." );

        static const double a1 = -3.969683028665376e+01;
        static const double a2 = 2.209460984245205e+02;
        static const double a3 = -2.759285104469687e+02;
        static const double a4 = 1.383577518672690e+02;
        static const double a5 = -3.066479806614716e+01;
        static const double a6 = 2.506628277459239e+00;

        static const double b1 = -5.447609879822406e+01;
        static const double b2 = 1.615858368580409e+02;
        static const double b3 = -1.556989798598866e+02;
        static const double b4 = 6.680131188771972e+01;
        static const double b5 = -1.328068155288572e+01;

        static const double c1 = -7.784894002430293e-03;
        static const double c2 = -3.223964580411365e-01;
        static const double c3 = -2.400758277161838e+00;
        static const double c4 = -2.549732539343734e+00;
        static const double c5 = 4.374664141464968e+00;
        static const double c6 = 2.938163982698783e+00;

        static const double d1 = 7.784695709041462e-03;
        static const double d2 = 3.224671290700398e-01;
        static const double d3 = 2.445134137142996e+00;
        static const double d4 = 3.754408661907416e+00;

        const double plow = 0.02425;
        const double phigh = 1 - plow;

        double q, r;
        if( p < plow )
        {
            // lower tail
            q = std::sqrt( -2 * std::log( p ) );
            return ( ( ( ( ( c1 * q + c2 ) * q + c3 ) * q + c4 ) * q + c5 ) * q
                       + c6 )
                   / ( ( ( ( d1 * q + d2 ) * q + d3 ) * q + d4 ) * q + 1 );
        }
        else if( p > phigh )
        {
            // upper tail
            q = std::sqrt( -2 * std::log( 1 - p ) );
            return -( ( ( ( ( c1 * q + c2 ) * q + c3 ) * q + c4 ) * q + c5 ) * q
                       + c6 )
                   / ( ( ( ( d1 * q + d2 ) * q + d3 ) * q + d4 ) * q + 1 );
        }

        // central region
        q = p - 0.5;
        r = q * q;
        return ( ( ( ( ( a1 * r + a2 ) * r + a3 ) * r + a4 ) * r + a5 ) * r
                   + a6 )
               * q
               / ( ( ( ( ( b1 * r + b2 ) * r + b3 ) * r + b4 ) * r + b5 ) * r
                   + 1 );
    }

    std::seed_seq create_seed_seq( uint64_t seed )
    {
        std::vector< uint32_t > seed_data = { static_cast< uint32_t >(
                                                  seed >> 32 ),
            static_cast< uint32_t >( seed & 0xFFFFFFFF ) };
        return std::seed_seq( seed_data.begin(), seed_data.end() );
    }

    // FNV-1a 64 bits hash: unlike absl::Hash, which is salted per process, it
    // gives the same seed for the same word across runs and platforms.
    uint64_t stable_hash( std::string_view word )
    {
        constexpr uint64_t FNV_OFFSET_BASIS{ 0xcbf29ce484222325 };
        constexpr uint64_t FNV_PRIME{ 0x100000001b3 };
        uint64_t hash{ FNV_OFFSET_BASIS };
        for( const auto character : word )
        {
            hash ^= static_cast< uint8_t >( character );
            hash *= FNV_PRIME;
        }
        return hash;
    }

    // Not named INFINITE: <windows.h> defines it as a macro
    constexpr auto DOUBLE_INFINITY = std::numeric_limits< double >::infinity();

    [[noreturn]] void throw_invalid_law( const std::string& law )
    {
        throw geode::OpenGeodeStochasticStochasticException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[RandomEngine] Cannot sample an invalid law: ", law };
    }

    // The rules of a law are in its is_valid(): they are not repeated here.
    // The exception is built in another function to keep this check inlined
    // in the sampling functions.
    template < typename Law >
    void check_law( const Law& law )
    {
        if( !law.is_valid() )
        {
            throw_invalid_law( law.string() );
        }
    }

    [[noreturn]] void throw_unrepresentable_sample( const std::string& law )
    {
        throw geode::OpenGeodeStochasticStochasticException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[RandomEngine] The law gives a sample out of the double range: ",
            law };
    }

    // A valid law can give values beyond the double range (huge mean or
    // standard deviation, exponent close to 1): raise an error instead of
    // returning inf or NaN.
    template < typename Law >
    double check_sample( double value, const Law& law )
    {
        if( !std::isfinite( value ) )
        {
            throw_unrepresentable_sample( law.string() );
        }
        return value;
    }

    uint64_t random_seed()
    {
        std::random_device device;
        return ( static_cast< uint64_t >( device() ) << 32 ) | device();
    }
} // namespace

namespace geode
{
    class RandomEngine::Impl
    {
    public:
        void set_seed( uint64_t number )
        {
            auto seed_seq = create_seed_seq( number );
            rand_gen_.seed( seed_seq );
        }

        void set_seed( std::string_view word )
        {
            set_seed( stable_hash( word ) );
        }

        template < typename Type >
        Type sample_uniform( const UniformClosed< Type >& law )
        {
            check_law( law );
            return absl::Uniform(
                absl::IntervalClosed, rand_gen_, law.min_value, law.max_value );
        }

        template < typename Type >
        Type sample_uniform( const UniformClosedOpen< Type >& law )
        {
            check_law( law );
            return absl::Uniform( absl::IntervalClosedOpen, rand_gen_,
                law.min_value, law.max_value );
        }

        double sample_gaussian( const Gaussian& law )
        {
            check_law( law );
            return check_sample(
                absl::gaussian_distribution< double >(
                    law.mean, law.standard_deviation )( rand_gen_ ),
                law );
        }

        double sample_truncated_gaussian( const TruncatedGaussian& law )
        {
            check_law( law );
            const auto min = law.min_value.value_or( -DOUBLE_INFINITY );
            const auto max = law.max_value.value_or( DOUBLE_INFINITY );
            const auto value =
                law.mean
                + law.standard_deviation
                      * sample_truncated_standard_normal(
                          ( min - law.mean ) / law.standard_deviation,
                          ( max - law.mean ) / law.standard_deviation );
            return check_sample( std::clamp( value, min, max ), law );
        }

        double sample_von_mises( const VonMises& law )
        {
            check_law( law );

            // Uniform approximation for very small concentration (nearly
            // uniform)
            if( law.concentration < geode::GLOBAL_EPSILON )
            {
                UniformClosedOpen< double > uniform_dist;
                uniform_dist.max_value = 2.0 * M_PI;
                return sample_uniform( uniform_dist );
            }

            // Normal approximation for very large concentration
            const double LARGE_KAPPA = 1e3; // threshold, can be tuned
            if( law.concentration > LARGE_KAPPA )
            {
                // Variance of approximate normal around mean
                const auto stddev = 1.0 / std::sqrt( law.concentration );
                absl::gaussian_distribution< double > normal_dist(
                    law.mean, stddev );
                return wrap_to_two_pi( normal_dist( rand_gen_ ) );
            }

            // Best & Fisher (1979) rejection algorithm for moderate
            // concentration
            const auto a =
                1.0
                + std::sqrt(
                    1.0 + 4.0 * law.concentration * law.concentration );
            const auto b =
                ( a - std::sqrt( 2.0 * a ) ) / ( 2.0 * law.concentration );
            const auto r = ( 1.0 + b * b ) / ( 2.0 * b );

            double theta;
            while( true )
            {
                auto u1 = unit_closed_( rand_gen_ );
                auto z = std::cos( M_PI * u1 );
                // z spans [-1, 1]: with |z|, f stayed positive and the
                // samples within +-pi/2 of the mean
                auto f = ( 1.0 + r * z ) / ( r + z );
                auto c = law.concentration * ( r - f );
                auto u2 = unit_closed_( rand_gen_ );

                if( u2 < c * ( 2.0 - c ) || u2 <= c * std::exp( 1.0 - c ) )
                {
                    theta = std::acos( f );
                    if( sample_bernoulli( 0.5 ) )
                        theta = -theta;
                    break;
                }
            }

            return wrap_to_two_pi( theta + law.mean );
        }

        double sample_axial_von_mises( const AxialVonMises& law )
        {
            check_law( law );
            VonMises doubled;
            doubled.mean = 2.0 * law.mean;
            doubled.concentration = law.concentration;
            // halving an angle of [0, 2pi[ gives an axis of [0, pi[
            return sample_von_mises( doubled ) / 2.0;
        }

        double sample_truncated_lognormal( const TruncatedLogNormal& law )
        {
            check_law( law );
            const auto min = law.min_value.value_or( 0. );
            const auto max = law.max_value.value_or( DOUBLE_INFINITY );
            // log( 0 ) = -inf and log( inf ) = inf are valid bounds
            const auto value =
                std::exp( law.mean
                          + law.standard_deviation
                                * sample_truncated_standard_normal(
                                    ( std::log( min ) - law.mean )
                                        / law.standard_deviation,
                                    ( std::log( max ) - law.mean )
                                        / law.standard_deviation ) );
            return check_sample( std::clamp( value, min, max ), law );
        }

        double sample_truncated_powerlaw( const TruncatedPowerLaw& law )
        {
            check_law( law );
            const auto xmin = law.lower_bound();
            const auto xmax = law.max_value.value_or( DOUBLE_INFINITY );
            absl::uniform_real_distribution< double > uniform( 0.0, 1.0 );
            const auto u = uniform( rand_gen_ );

            // Inverse CDF of the density proportional to x^-alpha on
            // [xmin, xmax], written with the ratio x / xmin so that no power
            // of xmin or xmax overflows. is_valid() ensures that xmax is
            // finite when alpha <= 1.
            if( std::fabs( law.alpha - 1.0 ) <= geode::GLOBAL_EPSILON
                && std::isfinite( xmax ) )
            {
                // bounded by the finite xmax: cannot overflow
                return std::clamp(
                    xmin * std::pow( xmax / xmin, u ), xmin, xmax );
            }
            const auto exponent = 1.0 - law.alpha;
            // with alpha > 1, an infinite xmax gives a null ratio
            const auto ratio = std::pow( xmax / xmin, exponent );
            return check_sample(
                std::clamp(
                    xmin * std::pow( 1.0 - u + u * ratio, 1.0 / exponent ),
                    xmin, xmax ),
                law );
        }

        double sample_log()
        {
            return std::log(
                absl::Uniform( absl::IntervalOpenClosed, rand_gen_, 0., 1.0 ) );
        }

        bool sample_bernoulli( double probability_of_success )
        {
            OpenGeodeStochasticStochasticException::check_exception(
                probability_of_success >= 0. && probability_of_success <= 1.0,
                nullptr, OpenGeodeException::TYPE::data,
                "[RandomEngine] Bernoulli sampling cannot be done since ",
                probability_of_success, " is not in [0.,1.]." );
            return absl::bernoulli_distribution( probability_of_success )(
                rand_gen_ );
        }

    private:
        // Samples the standard normal law truncated to [lower, upper] by
        // inverting its CDF. A range above the mean is mirrored below it,
        // where the CDF keeps its precision far in the tail.
        double sample_truncated_standard_normal( double lower, double upper )
        {
            if( lower > 0. )
            {
                return -sample_truncated_standard_normal( -upper, -lower );
            }
            // normal_quantile cannot handle 0 or 1
            const auto cdf_min =
                lower == -DOUBLE_INFINITY
                    ? std::numeric_limits< double >::min()
                    : std::max( normal_cdf( lower ),
                          std::numeric_limits< double >::min() );
            const auto cdf_max =
                upper == DOUBLE_INFINITY
                    ? 1. - std::numeric_limits< double >::epsilon()
                    : std::min( normal_cdf( upper ),
                          1. - std::numeric_limits< double >::epsilon() );
            OpenGeodeStochasticStochasticException::check_exception(
                cdf_min < cdf_max, nullptr, OpenGeodeException::TYPE::data,
                "[RandomEngine] Truncation range [", lower, ", ", upper,
                "] (in standard deviations) is too far in the tail of the "
                "normal law." );
            absl::uniform_real_distribution< double > uniform(
                cdf_min, cdf_max );
            return normal_quantile( uniform( rand_gen_ ) );
        }

        // std::mt19937_64 seeded by std::seed_seq: both are specified by the
        // C++ standard, so a seed gives the same sequence on every run and
        // platform (absl::BitGen salts its seed per process). Without seed,
        // the engine is randomly initialized. Values are drawn with the
        // absl:: distributions: the algorithms of the std:: ones are not
        // specified and differ between standard libraries.
        std::mt19937_64 rand_gen_{ random_seed() };
        // Uniform law on [0,1], built once: the same draws as sample_uniform
        // without validating the law nor computing nextafter( 1 ) each time.
        absl::uniform_real_distribution< double > unit_closed_{ 0.0,
            std::nextafter( 1.0, std::numeric_limits< double >::max() ) };
    };

    RandomEngine::RandomEngine() = default;

    RandomEngine::~RandomEngine() = default;

    void RandomEngine::set_seed( uint64_t number )
    {
        impl_->set_seed( number );
    }

    void RandomEngine::set_seed( std::string_view word )
    {
        impl_->set_seed( word );
    }

    template < typename Type >
    Type RandomEngine::sample_uniform( const UniformClosed< Type >& law )
    {
        return impl_->sample_uniform( law );
    }

    template index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosed< index_t >& );
    template local_index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosed< local_index_t >& );
    template signed_index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosed< signed_index_t >& );
    template float opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosed< float >& );
    template double opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosed< double >& );

    template < typename Type >
    Type RandomEngine::sample_uniform( const UniformClosedOpen< Type >& law )
    {
        return impl_->sample_uniform( law );
    }
    template index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosedOpen< index_t >& );
    template local_index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform(
            const UniformClosedOpen< local_index_t >& );
    template signed_index_t opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform(
            const UniformClosedOpen< signed_index_t >& );
    template float opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosedOpen< float >& );
    template double opengeode_stochastic_stochastic_api
        RandomEngine::sample_uniform( const UniformClosedOpen< double >& );

    double RandomEngine::sample_gaussian( const Gaussian& law )
    {
        return impl_->sample_gaussian( law );
    }
    double RandomEngine::sample_truncated_gaussian(
        const TruncatedGaussian& law )
    {
        return impl_->sample_truncated_gaussian( law );
    }
    double RandomEngine::sample_von_mises( const VonMises& law )
    {
        return impl_->sample_von_mises( law );
    }
    double RandomEngine::sample_axial_von_mises( const AxialVonMises& law )
    {
        return impl_->sample_axial_von_mises( law );
    }
    double RandomEngine::sample_truncated_lognormal(
        const TruncatedLogNormal& law )
    {
        return impl_->sample_truncated_lognormal( law );
    }
    double RandomEngine::sample_truncated_powerlaw(
        const TruncatedPowerLaw& law )
    {
        return impl_->sample_truncated_powerlaw( law );
    }
    double RandomEngine::sample_log()
    {
        return impl_->sample_log();
    }
    bool RandomEngine::sample_bernoulli( double probability_of_success )
    {
        return impl_->sample_bernoulli( probability_of_success );
    }

} // namespace geode