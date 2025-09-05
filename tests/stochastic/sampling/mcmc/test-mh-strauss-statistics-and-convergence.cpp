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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/sampling/direct/configuration_sampler/point_configuration_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/intensity_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
namespace
{

    void test_convergence( double domain_length,
        double poisson_density,
        double gamma,
        double nb_points,
        double nb_paires )
    {
        geode::Point2D min_point{ { 0., 0. } };
        geode::Point2D max_point{ { domain_length, domain_length } };

        geode::BoundingBox2D box;
        box.add_point( min_point );
        box.add_point( max_point );

        double area = domain_length * domain_length;

        geode::uuid group_id;
        geode::UniformPointConfigurationSampler< 2 > sampler( box, group_id );

        geode::GibbsEnergy< geode::Point2D > energy;
        energy.add_energy_term(
            std::make_unique< geode::IntensityTerm< geode::Point2D > >(
                poisson_density, group_id ) );

        auto interaction_fn = []( const geode::Point2D& a,
                                  const geode::Point2D& b ) {
            auto dx = a.value( 0 ) - b.value( 0 );
            auto dy = a.value( 1 ) - b.value( 1 );
            auto dist_sq = dx * dx + dy * dy;
            return dist_sq < 2;
        };

        geode::PairwiseTerm< geode::Point2D, decltype( interaction_fn ) > term(
            gamma, interaction_fn );
        energy.add_energy_term( std::make_unique<
            geode::PairwiseTerm< geode::Point2D, decltype( interaction_fn ) > >(
            gamma, interaction_fn ) );

        auto kernel1 =
            geode::create_birth_death_change_kernel< geode::Point2D >(
                sampler, 0.33, 0.33 );
        geode::MetropolisHastings< geode::Point2D > mh(
            energy, std::move( kernel1 ) );
        geode::RandomEngine engine;
        engine.set_seed( "@mh-test@" );
        std::unordered_map< geode::uuid, geode::index_t > targets = {
            { group_id, 0 }
        };
        geode::Configuration< geode::Point2D > state =
            mh.initialize_configuration_with_sampling( engine, targets );

        constexpr geode::index_t N{ 100000 };

        // Sampling
        double sum_points = 0.0;
        double sum_paires = 0.0;
        for( const auto i : geode::Range{ N } )
        {
            mh.step( state, engine );
            auto stats = energy.ordered_energy_term_statistics( state );
            sum_points += stats[0];
            sum_paires += stats[1];
        }

        double mean_points = sum_points / N;
        double mean_paires = sum_paires / N;

        geode::Logger::info( "[MH test] mean points = ", mean_points,
            " and mean paires = ", mean_paires, " (expected ", nb_points, " ; ",
            nb_paires, ") " );

        OPENGEODE_EXCEPTION( std::abs( mean_points - nb_points ) < 1,
            "[MH test] unexpected nb points." );

        // ------------------------------------------------------------
        // Variance test: Poisson => Var(N) ≈ E[N]
        // ------------------------------------------------------------
        if( nb_paires == 0 )
        {
            OPENGEODE_EXCEPTION(
                mean_paires == 0, "[MH test] unexpected nb paires." );
        }
        else
        {
            OPENGEODE_EXCEPTION( std::abs( mean_paires - nb_paires ) < 3,
                "[MH test] unexpected nb paires." );
        }
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        double domain_length{ 10. };
        double poisson_density_lambda{ 0.5 };
        std::array< double, 3 > gamma{ 0, 0.5, 1. };

        test_convergence( domain_length, poisson_density_lambda, 0., 15, 0 );
        test_convergence( domain_length, poisson_density_lambda, 1., 50, 70 );
        test_convergence( domain_length, poisson_density_lambda, 0.5, 24, 10 );

        geode::Logger::info( "MH STATISTICS AND CONVERGENCE TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}