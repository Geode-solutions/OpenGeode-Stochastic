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
#include <geode/stochastic/common.hpp>

#include <absl/random/bernoulli_distribution.h>
#include <absl/random/gaussian_distribution.h>
#include <absl/random/random.h>
#include <absl/random/uniform_int_distribution.h>
#include <absl/random/uniform_real_distribution.h>

#include <absl/hash/hash.h>
#include <limits>

namespace
{
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
        OPENGEODE_EXCEPTION( p >= 0.0 && p <= 1.0,
            "[normal_quantile] - p must be in (0,1). Check the consistencies "
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
        else
        {
            // central region
            q = p - 0.5;
            r = q * q;
            return ( ( ( ( ( a1 * r + a2 ) * r + a3 ) * r + a4 ) * r + a5 ) * r
                       + a6 )
                   * q
                   / ( ( ( ( ( b1 * r + b2 ) * r + b3 ) * r + b4 ) * r + b5 )
                           * r
                       + 1 );
        }
    }

    std::seed_seq create_seed_seq( uint64_t seed )
    {
        std::vector< uint32_t > seed_data = { static_cast< uint32_t >(
                                                  seed >> 32 ),
            static_cast< uint32_t >( seed & 0xFFFFFFFF ) };
        return std::seed_seq( seed_data.begin(), seed_data.end() );
    }

    std::seed_seq create_seed_seq( std::string_view word )
    {
        uint64_t seed = absl::Hash< std::string_view >{}( word );
        return create_seed_seq( seed );
    }
} // namespace

namespace geode
{
    class RandomEngine::Impl
    {
    public:
        void set_seed( uint64_t number )
        {
            rand_gen_ = absl::BitGen{ create_seed_seq( number ) };
        }

        void set_seed( std::string_view word )
        {
            rand_gen_ = absl::BitGen{ create_seed_seq( word ) };
        }

        template < typename Type >
        Type sample_uniform( const UniformClosed< Type >& law )
        {
            OPENGEODE_ASSERT( law.min_value <= law.max_value,
                "[Uniform sampling] - Wrong range ", law.min_value,
                " is not <= than ", law.max_value, "." );
            return absl::Uniform(
                absl::IntervalClosed, rand_gen_, law.min_value, law.max_value );
        }

        template < typename Type >
        Type sample_uniform( const UniformClosedOpen< Type >& law )
        {
            OPENGEODE_ASSERT( law.min_value < law.max_value,
                "[Uniform sampling] - Wrong range ", law.min_value,
                " is not < than ", law.max_value, "." );
            return absl::Uniform( absl::IntervalClosedOpen, rand_gen_,
                law.min_value, law.max_value );
        }

        double sample_gaussian( const Gaussian& law )
        {
            OPENGEODE_ASSERT( law.standard_deviation > 0
                                  && std::isfinite( law.standard_deviation )
                                  && std::isfinite( law.mean ),
                "[Gaussian sampling] - Infinite "
                "parameters or negative standard deviation N(",
                law.mean, law.standard_deviation, ")." );
            return absl::gaussian_distribution< double >(
                law.mean, law.standard_deviation )( rand_gen_ );
        }

        double sample_truncated_gaussian( const TruncatedGaussian& law )
        {
            OPENGEODE_ASSERT( law.standard_deviation > 0
                                  && std::isfinite( law.standard_deviation )
                                  && std::isfinite( law.mean ),
                "[Gaussian sampling] - Infinite parameters or "
                "negative standard deviation N(",
                law.mean, ",", law.standard_deviation, ")." );

            const double max = law.max_value.value_or(
                std::numeric_limits< double >::infinity() );
            const double min = law.min_value.value_or(
                -std::numeric_limits< double >::infinity() );

            OPENGEODE_ASSERT( min < max,
                "[Gaussian sampling] - Wrong truncation range ", min,
                " is not < than ", max, "." );

            // Standardize bounds
            const double alpha = ( min - law.mean ) / law.standard_deviation;
            const double beta = ( max - law.mean ) / law.standard_deviation;

            // Compute CDF of bounds, handling infinite alpha/beta
            double F_min = std::isfinite( alpha ) ? normal_cdf( alpha ) : 0.0;
            double F_max = std::isfinite( beta ) ? normal_cdf( beta ) : 1.0;

            // Clamp to avoid exact 0 or 1 (normal_quantile cannot handle them)
            F_min = std::max( F_min, geode::GLOBAL_EPSILON );
            F_max = std::min( F_max, 1.0 - geode::GLOBAL_EPSILON );

            OPENGEODE_EXCEPTION( F_min < F_max,
                "[Gaussian sampling] - truncation "
                "range is extreme please check inputs" );

            // Sample uniform in [F_min, F_max]
            std::uniform_real_distribution< double > uniform( F_min, F_max );
            const double u = uniform( rand_gen_ );

            // Map through inverse CDF
            return law.mean + law.standard_deviation * normal_quantile( u );
        }

        double sample_von_mises( const VonMises& law )
        {
            OPENGEODE_ASSERT(
                law.concentration >= 0.0 && std::isfinite( law.mean ),
                "[VonMises sampling] - Invalid parameters: mean=", law.mean,
                ", concentration=", law.concentration, "." );

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
                const double stddev = 1.0 / std::sqrt( law.concentration );
                std::normal_distribution< double > normal_dist(
                    law.mean, stddev );
                double theta = normal_dist( rand_gen_ );
                // Wrap to [0, 2π)
                return std::fmod( theta + 2.0 * M_PI, 2.0 * M_PI );
            }

            // Best & Fisher (1979) rejection algorithm for moderate
            // concentration
            const double a =
                1.0
                + std::sqrt(
                    1.0 + 4.0 * law.concentration * law.concentration );
            const double b =
                ( a - std::sqrt( 2.0 * a ) ) / ( 2.0 * law.concentration );
            const double r = ( 1.0 + b * b ) / ( 2.0 * b );

            double theta;
            UniformClosed< double > uniform_dist;
            while( true )
            {
                double u1 = sample_uniform( uniform_dist );
                double z = std::cos( M_PI * u1 );
                double f = ( 1.0 + r * std::abs( z ) ) / ( r + std::abs( z ) );
                double c = law.concentration * ( r - f );
                double u2 = sample_uniform( uniform_dist );

                if( u2 < c * ( 2.0 - c ) || u2 <= c * std::exp( 1.0 - c ) )
                {
                    theta = std::acos( f );
                    if( sample_bernoulli( 0.5 ) )
                        theta = -theta;
                    break;
                }
            }

            // Shift by mean and wrap to [0, 2π)
            theta += law.mean;
            theta = std::fmod( theta + 2.0 * M_PI, 2.0 * M_PI );

            return theta;
        }

        double sample_truncated_lognormal( const TruncatedLogNormal& law )
        {
            // Basic sanity checks
            OPENGEODE_ASSERT( law.standard_deviation > 0
                                  && std::isfinite( law.standard_deviation )
                                  && std::isfinite( law.mean ),
                "[Truncated LogNormal sampling] - Infinite parameters or "
                "negative standard deviation N(",
                law.mean, ", ", law.standard_deviation, ")." );

            // Determine min and max, respecting optional values
            const double min_val = law.min_value.value_or( 0.0 );
            const double max_val = law.max_value.value_or(
                std::numeric_limits< double >::infinity() );

            OPENGEODE_ASSERT( min_val < max_val,
                "[Truncated LogNormal sampling] - Wrong truncation range ",
                min_val, " is not < than ", max_val, "." );

            // Transform to standard normal space
            const double zmin =
                ( std::log( min_val ) - law.mean ) / law.standard_deviation;
            const double zmax =
                ( std::isfinite( max_val )
                        ? ( std::log( max_val ) - law.mean )
                              / law.standard_deviation
                        : std::numeric_limits< double >::infinity() );

            // Compute CDF bounds, handling infinite zmin/zmax
            double F_min =
                std::max( normal_cdf( zmin ), geode::GLOBAL_EPSILON );
            double F_max =
                std::min( normal_cdf( zmax ), 1.0 - geode::GLOBAL_EPSILON );

            OPENGEODE_EXCEPTION( F_min < F_max,
                "[Truncated LogNormal sampling] - truncation "
                "range is extreme please chack inputs" );

            // Sample uniform in [Fmin, Fmax]
            std::uniform_real_distribution< double > uniform( F_min, F_max );
            const double u = uniform( rand_gen_ );

            // Map through inverse CDF and exponentiate
            const double z = normal_quantile( u );
            return std::exp( law.mean + law.standard_deviation * z );
        }

        double sample_truncated_powerlaw( const TruncatedPowerLaw& law )
        {
            OPENGEODE_ASSERT(
                law.alpha > 0, "Power-law exponent alpha must be > 0" );

            // Set bounds
            const double xmin = law.min_value.value_or(
                geode::GLOBAL_EPSILON ); // default 1.0 if unspecified
            const double xmax = law.max_value.value_or(
                std::numeric_limits< double >::infinity() );

            OPENGEODE_ASSERT( xmin < xmax, "Truncated power-law: min >= max" );

            // Sample uniform
            std::uniform_real_distribution< double > uniform( 0.0, 1.0 );
            const double u = uniform( rand_gen_ );

            // Inverse CDF
            if( std::abs( law.alpha - 1.0 ) > geode::GLOBAL_EPSILON )
            {
                double xmin_pow = std::pow( xmin, 1.0 - law.alpha );
                double xmax_pow =
                    std::isfinite( xmax )
                        ? std::pow( xmax, 1.0 - law.alpha )
                        : std::numeric_limits< double >::infinity();
                if( !std::isfinite( xmax_pow ) )
                {
                    return std::pow( xmin_pow + u, 1.0 / ( 1.0 - law.alpha ) );
                }
                return std::pow( u * ( xmax_pow - xmin_pow ) + xmin_pow,
                    1.0 / ( 1.0 - law.alpha ) );
            }
            else
            {
                // alpha == 1
                const double xmax_eff =
                    std::isfinite( xmax ) ? xmax : xmin * geode::GLOBAL_EPSILON;
                return xmin * std::pow( xmax_eff / xmin, u );
            }
        }

        double sample_log()
        {
            return std::log(
                absl::Uniform( absl::IntervalOpenClosed, rand_gen_, 0., 1.0 ) );
        }

        bool sample_bernoulli( double probability_of_success )
        {
            OPENGEODE_ASSERT(
                probability_of_success >= 0. && probability_of_success <= 1.0,
                "Bernoulli sampling cannot be done since ",
                probability_of_success, " is not in [0.,1.]." );
            return absl::bernoulli_distribution( probability_of_success )(
                rand_gen_ );
        }

    private:
        absl::BitGen rand_gen_;
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