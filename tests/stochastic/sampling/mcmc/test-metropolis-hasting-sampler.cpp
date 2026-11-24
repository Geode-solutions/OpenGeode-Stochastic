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
#include <cmath>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/applications/poisson_process.hpp>
#include <geode/stochastic/inference/statistics_tools.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace
{
    constexpr double DOMAIN_SIZE{ 10. };

    geode::PoissonProcessBuilder< geode::Point2D > poisson_builder(
        double lambda, std::optional< double > expected_count )
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        domain.buffer_size = 0.;
        geode::PoissonProcessBuilder< geode::Point2D > poisson;
        poisson.set_domain( domain );
        auto& set = poisson.add_set( "points", lambda, expected_count );
        // NOLINTBEGIN(*-magic-numbers)
        set.fixed_objects = { geode::Point2D{ { 1., 1. } },
            geode::Point2D{ { 5., 5. } } };
        // NOLINTEND(*-magic-numbers)
        return poisson;
    }

    void test_beta_validation()
    {
        geode::Logger::info( "TEST - MH BETA VALIDATION" );
        // NOLINTBEGIN(*-magic-numbers)
        auto context =
            poisson_builder( 0.5, std::nullopt ).build_simulation_context();
        auto& mh = *context.mh_sampler;
        geode::OpenGeodeStochasticStochasticException::test(
            mh.beta() == 1., "[MH test] The default beta should be 1." );
        mh.set_beta( 0.5 );
        geode::OpenGeodeStochasticStochasticException::test(
            mh.beta() == 0.5, "[MH test] beta not set." );
        bool exception_thrown{ false };
        try
        {
            mh.set_beta( -1. );
        }
        catch( const geode::OpenGeodeException& )
        {
            exception_thrown = true;
        }
        geode::OpenGeodeStochasticStochasticException::test(
            exception_thrown && mh.beta() == 0.5,
            "[MH test] A negative beta should be refused." );
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    // Each step result must match the change of the state: an accepted birth
    // adds a free object, an accepted death removes one, a change or a
    // rejection keeps the count, and the fixed objects are never touched.
    void test_step_bookkeeping()
    {
        geode::Logger::info( "TEST - MH STEP BOOKKEEPING" );
        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-step-bookkeeping@" );
        // NOLINTBEGIN(*-magic-numbers)
        auto context =
            poisson_builder( 0.5, std::nullopt ).build_simulation_context();
        auto& state = *context.object_sets;
        const auto& const_state = state;
        const auto& set =
            const_state.get_set( const_state.get_set_uuid( "points" ) );
        const auto nb_fixed = set.nb_fixed_objects();
        geode::index_t accepted[3]{ 0, 0, 0 };
        for( const auto step : geode::Range{ 20000 } )
        {
            geode_unused( step );
            const auto nb_free_before = set.nb_free_objects();
            const auto result = context.mh_sampler->step( state, engine );
            const auto nb_free_after = set.nb_free_objects();
            geode::OpenGeodeStochasticStochasticException::test(
                set.nb_fixed_objects() == nb_fixed,
                "[MH test] The fixed objects should not change." );
            if( result.move_type == geode::MOVE_TYPE::invalid )
            {
                // death proposed on a set without free object: nothing is
                // decided and the state is unchanged
                geode::OpenGeodeStochasticStochasticException::test(
                    result.decision == geode::MH_DECISION::undecided
                        && nb_free_before == 0 && nb_free_after == 0,
                    "[MH test] An invalid move should only happen without "
                    "free object, and change nothing." );
                continue;
            }
            geode::OpenGeodeStochasticStochasticException::test(
                result.decision != geode::MH_DECISION::undecided,
                "[MH test] A valid move should be decided." );
            geode::OpenGeodeStochasticStochasticException::test(
                !std::isnan( result.log_accept ),
                "[MH test] The acceptance should be a number." );
            auto expected = nb_free_before;
            if( result.decision == geode::MH_DECISION::accepted )
            {
                if( result.move_type == geode::MOVE_TYPE::birth )
                {
                    expected++;
                    accepted[0]++;
                }
                else if( result.move_type == geode::MOVE_TYPE::death )
                {
                    expected--;
                    accepted[1]++;
                }
                else
                {
                    accepted[2]++;
                }
            }
            geode::OpenGeodeStochasticStochasticException::test(
                nb_free_after == expected,
                "[MH test] The number of free objects does not match the "
                "step result." );
        }
        geode::Logger::info( "Accepted moves: birth=", accepted[0],
            " death=", accepted[1], " change=", accepted[2] );
        geode::OpenGeodeStochasticStochasticException::test(
            accepted[0] > 0 && accepted[1] > 0 && accepted[2] > 0,
            "[MH test] Every move type should have been accepted." );
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    // The tempered chain targets a density proportional to exp(-beta * U):
    // for a Poisson process of intensity lambda, it is a Poisson process of
    // intensity lambda^beta.
    void test_beta_tempered_poisson()
    {
        geode::Logger::info( "TEST - MH TEMPERED POISSON" );
        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-tempered-poisson@" );
        // NOLINTBEGIN(*-magic-numbers)
        constexpr double lambda{ 0.25 };
        constexpr double nb_fixed{ 2. };
        for( const double beta : { 0., 0.5, 1. } )
        {
            const auto expected_count =
                std::pow( lambda, beta ) * DOMAIN_SIZE * DOMAIN_SIZE + nb_fixed;
            auto poisson = poisson_builder( lambda, expected_count );
            auto context = poisson.build_simulation_context();
            context.mh_sampler->set_beta( beta );
            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                context ) };
            geode::SimulationConfigurator config;
            config.realizations = 2000;
            config.metropolis_hasting_steps = 100;
            config.burn_in_steps = 2000;
            const auto tracker = runner.run( engine, config );
            const geode::TargetStatistics targets{ runner.model(),
                poisson.expected_statistics() };
            for( const auto& term_uuid : targets.active_terms() )
            {
                geode::Logger::info( "beta=", beta,
                    " mean statistic=", tracker.mean( term_uuid ),
                    " (target: ", targets.target( term_uuid ), ")" );
            }
            geode::statistics::validate( tracker, targets );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_beta_validation();
        test_step_bookkeeping();
        test_beta_tempered_poisson();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
