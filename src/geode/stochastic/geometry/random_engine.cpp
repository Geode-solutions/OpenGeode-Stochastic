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

#include <geode/stochastic/geometry/random_engine.hpp>

#include <geode/stochastic/geometry/distributions.hpp>

#include <geode/basic/pimpl_impl.hpp>

#include <absl/random/bernoulli_distribution.h>
#include <absl/random/gaussian_distribution.h>
#include <absl/random/random.h>
#include <absl/random/uniform_int_distribution.h>
#include <absl/random/uniform_real_distribution.h>

#include <absl/hash/hash.h>
#include <limits>
#include <variant>

namespace
{
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

        double sample( const Distribution& dist )
        {
            return std::visit(
                [this]( auto&& d ) {
                    return sample_distribution( d );
                },
                dist );
        }
        double sample_distribution( const Uniform< double >& law )
        {
            return absl::Uniform( absl::IntervalClosedOpen, rand_gen_,
                law.min.value, law.max.value );
        }
        double sample_distribution( const Gaussian& law )
        {
            return sample_gaussian( law );
        }
        template < typename Type >
        Type sample_uniform( const Uniform< Type >& law )
        {
            OPENGEODE_DATA_EXCEPTION( law.min.value < law.max.value,
                "Uniform sampling cannot be done since ", law.min.value,
                " is not < than ", law.max.value, "." );
            if( law.min.is_included )
            {
                if( law.max.is_included )
                {
                    return absl::Uniform( absl::IntervalClosed, rand_gen_,
                        law.min.value, law.max.value );
                }
                return absl::Uniform( absl::IntervalClosedOpen, rand_gen_,
                    law.min.value, law.max.value );
            }
            if( law.max.is_included )
            {
                return absl::Uniform( absl::IntervalOpenClosed, rand_gen_,
                    law.min.value, law.max.value );
            }
            return absl::Uniform(
                absl::IntervalOpen, rand_gen_, law.min.value, law.max.value );
        }

        double sample_gaussian( const Gaussian& law )
        {
            OPENGEODE_DATA_EXCEPTION(
                law.standard_deviation > 0
                    && std::isfinite( law.standard_deviation )
                    && std::isfinite( law.mean ),
                "Gaussian sampling cannot be done, please check the mean (",
                law.mean, ") and the standard deviation (",
                law.standard_deviation, ")." );
            if( !law.min.has_value() || !law.max.has_value() )
            {
                return absl::gaussian_distribution< double >(
                    law.mean, law.standard_deviation )( rand_gen_ );
            };
            const auto max =
                law.max.value_or( std::numeric_limits< double >::infinity() );
            const auto min =
                law.min.value_or( -std::numeric_limits< double >::infinity() );

            OPENGEODE_DATA_EXCEPTION( min < max,
                "Gaussian sampling cannot be done since ", min,
                " is not < than ", max, "." );

            OPENGEODE_DATA_EXCEPTION( law.mean > min && law.mean < max,
                "Gaussian sampling cannot be done since please the mean value "
                "(",
                law.mean, ") is not in [", min, ",", max, "]." );
            double value;
            do
            {
                value = absl::gaussian_distribution< double >(
                    law.mean, law.standard_deviation )( rand_gen_ );
            } while( value < min || value > max );
            return value;
        }

        bool sample_bernoulli( double probability_of_success )
        {
            OPENGEODE_DATA_EXCEPTION(
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

    double RandomEngine::sample( const Distribution& dist )
    {
        return impl_->sample( dist );
    }

    template < typename Type >
    Type RandomEngine::sample_uniform( const Uniform< Type >& law )
    {
        return impl_->sample_uniform( law );
    }

    template index_t opengeode_stochastic_geometry_api
        RandomEngine::sample_uniform( const Uniform< index_t >& );
    template local_index_t opengeode_stochastic_geometry_api
        RandomEngine::sample_uniform( const Uniform< local_index_t >& );
    template signed_index_t opengeode_stochastic_geometry_api
        RandomEngine::sample_uniform( const Uniform< signed_index_t >& );
    template float opengeode_stochastic_geometry_api
        RandomEngine::sample_uniform( const Uniform< float >& );
    template double opengeode_stochastic_geometry_api
        RandomEngine::sample_uniform( const Uniform< double >& );

    double RandomEngine::sample_gaussian( const Gaussian& law )
    {
        return impl_->sample_gaussian( law );
    }

    bool RandomEngine::sample_bernoulli( double probability_of_success )
    {
        return impl_->sample_bernoulli( probability_of_success );
    }

} // namespace geode