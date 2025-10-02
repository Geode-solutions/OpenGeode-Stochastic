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
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/density_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/spatial/object_set.hpp>
namespace
{

    // ------------------------------------------------------------
    // Convergence test: mean number of points ≈ λ × area
    // ------------------------------------------------------------
    void test_convergence_one_subsete(
        geode::MetropolisHastings< geode::Point2D >& mh,
        const absl::flat_hash_map< geode::uuid, geode::index_t >& group_targets,
        geode::RandomEngine& engine,
        const geode::GibbsEnergy< geode::Point2D >& energy )
    {
        geode::ObjectSet< geode::Point2D > state =
            mh.initialize_object_set_with_sampling( engine, group_targets );
        mh.walk( state, engine, 500 );
        constexpr geode::index_t N{ 5000 };

        // Sampling

        double sum_points = 0.0;
        double sum_sq = 0.0;
        for( const auto i : geode::Range{ N } )
        {
            mh.walk( state, engine, 1000 );
            auto stats = energy.ordered_energy_term_statistics( state );
            sum_points += stats[0];
            sum_sq += stats[0] * stats[0];
        }

        double mean_points = sum_points / N;
        double var = ( sum_sq / N ) - ( mean_points * mean_points );

        for( const auto& [group_uuid, expected_points] : group_targets )
        {
            geode_unused( group_uuid );
            geode::Logger::info( "[MH test] mean points = ", mean_points,
                " and var = ", var, " (expected ", expected_points, ")" );

            OPENGEODE_EXCEPTION( std::abs( mean_points - expected_points )
                                     < 0.01 * expected_points,
                "[MH test] mean number of points not close to expected." );

            // ------------------------------------------------------------
            // Variance test: Poisson => Var(N) ≈ E[N]
            // ------------------------------------------------------------
            OPENGEODE_EXCEPTION(
                std::abs( var - mean_points ) < 0.1 * expected_points,
                "[MH test] variance not close to Poisson expectation." );
        }
    }

    // ------------------------------------------------------------
    // Dynamic independence: different kernels → same stationary law
    // ------------------------------------------------------------
    void test_dynamic_independence(
        double domain_length, double poisson_density, double birth_ratio )
    {
        geode::RandomEngine engine;
        engine.set_seed( "@mh-test@" );

        geode::Point2D min_point{ { 0., 0. } };
        geode::Point2D max_point{ { domain_length, domain_length } };

        geode::BoundingBox2D box;
        box.add_point( min_point );
        box.add_point( max_point );

        double area = domain_length * domain_length;
        geode::uuid subset_id;
        geode::UniformPointSetSampler< 2 > sampler( box, subset_id );

        geode::GibbsEnergy< geode::Point2D > poisson_energy;
        poisson_energy.add_energy_term(
            std::make_unique< geode::DensityTerm< geode::Point2D > >(
                "intensity", poisson_density, subset_id ) );

        absl::flat_hash_map< geode::uuid, geode::index_t > targets = {
            { subset_id,
                static_cast< geode::index_t >( poisson_density * area ) }
        };
        // Kernel with only birth/death
        auto kernel1 = geode::create_birth_death_kernel< geode::Point2D >(
            sampler, birth_ratio );
        geode::MetropolisHastings< geode::Point2D > mh1(
            poisson_energy, std::move( kernel1 ) );

        // Kernel with birth/death/change
        auto kernel2 =
            geode::create_birth_death_change_kernel< geode::Point2D >(
                sampler, birth_ratio, 0.1 );
        geode::MetropolisHastings< geode::Point2D > mh2(
            poisson_energy, std::move( kernel2 ) );

        geode::Logger::info(
            "[MH test] Testing kernel with birth/death only..." );
        test_convergence_one_subsete( mh1, targets, engine, poisson_energy );

        geode::Logger::info(
            "[MH test] Testing kernel with birth/death/change..." );
        test_convergence_one_subsete( mh2, targets, engine, poisson_energy );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        double domain_length{ 10. };
        double poisson_density_lambda{ 0.5 };
        std::array< double, 7 > birth_ratio{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
            0.8 };
        for( const auto b : birth_ratio )
        {
            test_dynamic_independence(
                domain_length, poisson_density_lambda, b );
        }

        geode::Logger::info( "MH STATISTICS AND CONVERGENCE TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}