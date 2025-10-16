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
#include <geode/stochastic/spatial/object_sets.hpp>
namespace
{
    struct PoissonDescription
    {
        std::string name;
        double density;

        // mh dynamic
        double death_birth_ratio{ 2. };
        double birth_ratio{ 0.5 };
        double change_ratio{ 1. };
    };

    struct MultitypePoissonDescription
    {
        // voi
        geode::Point2D min_point;
        geode::Point2D max_point;

        // object sets
        std::vector< PoissonDescription > set_desc;

        // mh
        geode::index_t nb_steps{ 1000 };
        geode::index_t nb_realizations{ 1000 };
    };
    struct UserProblem
    {
        geode::BoundingBox2D box;
        geode::ObjectSets< geode::Point2D > object_set;
        std::vector< geode::uuid > object_set_id;

        geode::GibbsEnergy< geode::Point2D > gibbs_energy;

        absl::flat_hash_map< geode::uuid, geode::uuid > set_energy_term_ids;

        absl::flat_hash_map< geode::uuid, double > set_stats_targets;

        absl::flat_hash_map< geode::uuid, geode::UniformPointSetSampler< 2 > >
            set_samplers;

        std::unique_ptr< geode::MetropolisHastings< geode::Point2D > >
            mh_sampler;
    };

    UserProblem create_problems(
        const MultitypePoissonDescription& description )
    {
        UserProblem problem;
        problem.box.add_point( description.min_point );
        problem.box.add_point( description.max_point );
        double area = problem.box.n_volume();

        std::unique_ptr< geode::ProposalKernel< geode::Point2D > >
            proposal_kernel =
                std::make_unique< geode::ProposalKernel< geode::Point2D > >();
        for( const auto& points_desc : description.set_desc )
        {
            auto set_id = problem.object_set.add_set();
            problem.object_set_id.push_back( set_id );

            // this should be linked to the object subset
            // flat_hash_map<set_id,ObjectSetSampler>
            problem.set_samplers.emplace( set_id,
                geode::UniformPointSetSampler< 2 >{ problem.box, set_id } );
            OPENGEODE_EXCEPTION( points_desc.death_birth_ratio > 0.,
                "Object cannot be add or removed. Please set a BIRTH-DEATH "
                "with a positive probability." );
            proposal_kernel->add_move(
                std::make_unique< geode::BirthDeathMove< geode::Point2D > >(
                    problem.set_samplers.at( set_id ),
                    points_desc.death_birth_ratio, points_desc.birth_ratio ) );
            if( points_desc.change_ratio > 0. )
            {
                proposal_kernel->add_move(
                    std::make_unique< geode::ChangeMove< geode::Point2D > >(
                        problem.set_samplers.at( set_id ),
                        points_desc.change_ratio ) );
            }

            // energy terms here we define intra subset terms (can be
            // several of them) need to add inter set energy terms
            problem.set_energy_term_ids.emplace( set_id,
                problem.gibbs_energy.add_energy_term(
                    std::make_unique< geode::DensityTerm< geode::Point2D > >(
                        points_desc.name, points_desc.density, set_id ) ) );

            problem.set_stats_targets.emplace(
                set_id, points_desc.density * area );
        }
        problem.mh_sampler =
            std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                problem.gibbs_energy, std::move( proposal_kernel ) );

        return problem;
    }
    void test_convergence(
        const MultitypePoissonDescription& problem_description,
        geode::RandomEngine& engine )
    {
        auto problem = create_problems( problem_description );

        problem.mh_sampler->walk( problem.object_set, engine, 500 );

        std::vector< double > sum_points( problem.object_set_id.size(), 0. );
        std::vector< double > sum_sq_points( problem.object_set_id.size(), 0. );
        auto N = problem_description.nb_realizations;

        for( const auto i : geode::Range{ N } )
        {
            problem.mh_sampler->walk(
                problem.object_set, engine, problem_description.nb_steps );
            for( const auto set_id :
                geode::Range{ problem.object_set_id.size() } )
            {
                const auto& energy_term_uuid = problem.set_energy_term_ids.at(
                    problem.object_set_id[set_id] );
                const auto nb_points =
                    problem.gibbs_energy.energy_term_statistic(
                        problem.object_set, energy_term_uuid );
                sum_points[set_id] += nb_points;
                sum_sq_points[set_id] += nb_points * nb_points;
            }
        }
        std::transform( sum_points.begin(), sum_points.end(),
            sum_points.begin(), [N]( double p ) {
                return p / N;
            } );
        std::transform( sum_sq_points.begin(), sum_sq_points.end(),
            sum_sq_points.begin(), [N]( double p ) {
                return p / N;
            } );
        for( const auto set_id : geode::Range{ problem.object_set_id.size() } )
        {
            const auto& set_id = problem.object_set_id[set_id];
            const auto expected_points = problem.set_stats_targets.at( set_id );

            const auto variance = sum_sq_points[set_id]
                                  - ( sum_points[set_id] * sum_points[set_id] );
            geode::Logger::info( "[MH test] mean points = ", sum_points[set_id],
                " and var = ", variance, " (expected ", expected_points, ")" );
            const auto error_stat =
                std::abs( sum_points[set_id] - expected_points )
                / expected_points;
            OPENGEODE_EXCEPTION( error_stat < 0.02,
                "[MH test] mean number of points not close to enought to "
                "expected value --> error : ",
                error_stat );
            // Variance test: Poisson => Var(N) ≈ E[N]
            const auto error_var =
                std::abs( variance - expected_points ) / expected_points;
            OPENGEODE_EXCEPTION( error_var < 0.15,
                "[MH test] Variance not close to enought to "
                "expected value --> error : ",
                error_var );
        }
    }

    void test_single_type_poisson()
    {
        geode::Logger::info( "TEST - MH SINGLE TYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-POISSON@" );

        std::array< double, 4 > birth_ratio{ 0.1, 0.3, 0.7, 0.9 };
        std::array< double, 4 > change_ratio{ 0., 1., 1., 0. };

        for( const auto config : geode::Range{ birth_ratio.size() } )
        {
            MultitypePoissonDescription problem_description;
            problem_description.min_point = geode::Point2D{ { 0., 0. } };
            problem_description.max_point = geode::Point2D{ { 10., 10. } };

            std::vector< PoissonDescription > poisson_description{ { "set01",
                0.5, 1., birth_ratio[config], change_ratio[config] } };
            problem_description.set_desc = poisson_description;

            problem_description.nb_steps = 1000.;
            problem_description.nb_realizations = 1000.;

            test_convergence( problem_description, engine );
        }
        geode::Logger::info( "MH SINGLE TYPE POISSON -- SUCCESS!" );
    }

    void test_multitype_poisson()
    {
        geode::Logger::info( "TEST - MH MULTITYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-POISSON-multi@" );

        MultitypePoissonDescription problem_description;
        problem_description.min_point = geode::Point2D{ { 0., 0. } };
        problem_description.max_point = geode::Point2D{ { 10., 10. } };

        std::vector< PoissonDescription > poisson_description{
            { "set01", 0.1, 3., 0.25, 1. }, { "set02", 0.4, 1., 0.75, 1. },
            { "set03", 0.3, 5., 0.5, 1. }
        };
        problem_description.set_desc = poisson_description;

        problem_description.nb_steps = 1000.;
        problem_description.nb_realizations = 1000.;

        test_convergence( problem_description, engine );
        geode::Logger::info( "MH MULTITYPE POISSON -- SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        test_single_type_poisson();
        test_multitype_poisson();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}