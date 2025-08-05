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

#include <geode/stochastic/sampling/direct/ball_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/distance.hpp>

const int NUMBER_OF_STEPS = 10000;

// ------------------------- MCMCMove Interface -------------------------
template < typename State >
class MCMCMove
{
public:
    virtual ~MCMCMove() = default;
    virtual State propose(
        const State& current, geode::RandomEngine& rng ) const = 0;
    virtual double acceptance_log_ratio( const State&, const State& ) const = 0;
};

// ------------------------- Gaussian Perturbation Move
// -------------------------
class GaussianMove : public MCMCMove< std::vector< double > >
{
public:
    explicit GaussianMove( double stddev ) : stddev_( stddev ) {}

    std::vector< double > propose( const std::vector< double >& current,
        geode::RandomEngine& rng ) const override
    {
        std::vector< double > next = current;
        for( auto& x : next )
        {
            geode::Gaussian dist{ 0.0, stddev_ };
            x += rng.sample_gaussian( dist );
        }
        return next;
    }

    double acceptance_log_ratio( const std::vector< double >&,
        const std::vector< double >& ) const override
    {
        return 0.0; // symmetric
    }

private:
    double stddev_;
};
class MCMCSampler
{
public:
    MCMCSampler( std::vector< double > initial,
        std::function< double( const std::vector< double >& ) > log_prob,
        geode::RandomEngine& rng )
        : state_( std::move( initial ) ),
          log_prob_( std::move( log_prob ) ),
          rng_( rng )
    {
    }

    void set_move( std::unique_ptr< MCMCMove< std::vector< double > > > move )
    {
        move_ = std::move( move );
    }

    void step()
    {
        auto proposed = move_->propose( state_, rng_ );
        double log_accept_ratio =
            log_prob_( proposed ) - log_prob_( state_ )
            + move_->acceptance_log_ratio( state_, proposed );

        if( rng_.sample_bernoulli(
                std::min( 1.0, std::exp( log_accept_ratio ) ) ) )
        {
            state_ = std::move( proposed );
        }
    }

    const std::vector< double >& current() const
    {
        return state_;
    }

private:
    std::vector< double > state_;
    std::function< double( const std::vector< double >& ) > log_prob_;
    std::unique_ptr< MCMCMove< std::vector< double > > > move_;
    geode::RandomEngine& rng_;
};

double log_target( const std::vector< double >& x )
{
    double sum = 0.0;
    for( double xi : x )
    {
        sum -= 0.5 * xi * xi; // standard normal log density (ignoring const)
    }
    return sum;
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        constexpr int dim = 2;
        std::vector< double > initial( dim, 0.0 );
        geode::RandomEngine rng;
        rng.set_seed( 42 );

        MCMCSampler sampler( initial, log_target, rng );
        sampler.set_move( std::make_unique< GaussianMove >( 0.5 ) );

        for( int i = 0; i < NUMBER_OF_STEPS; ++i )
        {
            sampler.step();
            const auto& state = sampler.current();
            auto str = absl::StrCat( "Step ", i, ": (" );
            for( auto v : state )
                absl::StrAppend( &str, v, " " );
            geode::Logger::info( str, ")" );
        }
        geode::Logger::info( "TEST SUCCESS" );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
