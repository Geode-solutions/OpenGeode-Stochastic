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

#include <geode/stochastic/geometry/random_sample_spec.hpp>

#include <geode/basic/pimpl_impl.hpp>

#include "absl/random/bernoulli_distribution.h"
#include "absl/random/gaussian_distribution.h"
#include "absl/random/random.h"
#include "absl/random/uniform_int_distribution.h"
#include "absl/random/uniform_real_distribution.h"

#include <variant>

namespace
{
    std::seed_seq create_seed_seq( uint64_t seed )
    {
        std::vector< uint32_t > seed_data = { static_cast< uint32_t >(
                                                  seed & 0xFFFFFFFF ),
            static_cast< uint32_t >( seed >> 32 ) };
        return std::seed_seq( seed_data.begin(), seed_data.end() );
    }
} // namespace

namespace geode
{
    class RandomEngine::Impl
    {
    public:
        explicit Impl( uint64_t seed ) : rand_gen_{ create_seed_seq( seed ) } {}
        ~Impl() = default;

        template < typename RANDSPEC >
        auto do_sample( const RANDSPEC& spec )
        {
            return do_sample_impl( spec );
        }

    private:
        int do_sample_impl( const UniformInt& spec )
        {
            return absl::Uniform(
                absl::IntervalClosed, rand_gen_, spec.min, spec.max );
        }
        double do_sample_impl( const UniformDouble& spec )
        {
            return absl::Uniform( rand_gen_, spec.min, spec.max );
        }
        double do_sample_impl( const Gaussian& spec )
        {
            return absl::gaussian_distribution< double >(
                spec.mean, spec.stddev )( rand_gen_ );
        }
        bool do_sample_impl( const Bernoulli& spec )
        {
            return absl::bernoulli_distribution( spec.probability )(
                rand_gen_ );
        }

    private:
        absl::BitGen rand_gen_;
    };

    RandomEngine::RandomEngine( uint64_t seed ) : impl_{ seed } {}
    RandomEngine::~RandomEngine() = default;

    double RandomEngine::sample( const RandomDoubleSpec& random_double_spec )
    {
        return std::visit(
            [this]( const auto& s ) {
                return impl_->do_sample( s );
            },
            random_double_spec );
    }

    template < typename RANDSPEC >
    auto RandomEngine::sample( const RANDSPEC& spec )
    {
        return impl_->do_sample( spec );
    }

    template auto RandomEngine::sample( const UniformInt& );
    template auto RandomEngine::sample( const UniformDouble& );
    template auto RandomEngine::sample( const Gaussian& );
    template auto RandomEngine::sample( const Bernoulli& );

} // namespace geode