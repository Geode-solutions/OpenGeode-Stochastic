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
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/spatial/object_set.hpp>

namespace
{
    struct StraussDescription
    {
        std::string name;
        double density;
        double expected_number_of_objects;

        // interaction
        geode::PairwiseInteraction< geode::Point2D >::SCOPE scope;
        double distance_treshold;
        double gamma;
        double expected_number_of_intersections;

        // mh dynamic
        double death_birth_ratio{ 1. };
        double birth_ratio{ 0.5 };
        double change_ratio{ 1. };
    };

    struct MultitypeStraussDescription
    {
        // voi
        geode::Point2D min_point;
        geode::Point2D max_point;

        // object sets
        std::vector< StraussDescription > set_desc;

        // mh
        geode::index_t nb_steps{ 10000 };
        geode::index_t nb_realizations{ 1000 };
    };
    struct UserProblem
    {
        geode::BoundingBox2D box;
        geode::ObjectSet< geode::Point2D > object_set;
        std::vector< geode::uuid > object_subset_id;

        geode::GibbsEnergy< geode::Point2D > gibbs_energy;

        absl::flat_hash_map< geode::uuid, std::vector< geode::uuid > >
            subset_energy_term_ids;
        absl::flat_hash_map< geode::uuid, std::vector< double > >
            subset_stats_targets;

        absl::flat_hash_map< geode::uuid, geode::UniformPointSetSampler< 2 > >
            subset_samplers;
        std::unique_ptr< geode::MetropolisHastings< geode::Point2D > >
            mh_sampler;

        std::string string() const
        {
            std::string message( "User Problem for stochastic simulation: " );
            absl::StrAppend(
                &message, "\n\t - VOI: BoundingBox ", box.string() );
            absl::StrAppend( &message, "\n\t - ", object_set.string() );
            absl::StrAppend( &message, "\n\t - subset uuid list: " );
            for( const auto& uuid : object_subset_id )
            {
                absl::StrAppend(
                    &message, "\n\t --> subset uuid: ", uuid.string() );
            }
            absl::StrAppend( &message, "\n\t - ", gibbs_energy.string() );
            absl::StrAppend(
                &message, "\n\t - subset energy term uuid list: " );
            for( const auto& [uuid, et_uuids] : subset_energy_term_ids )
            {
                absl::StrAppend(
                    &message, "\n\t --> subset uuid: ", uuid.string() );
                for( const auto& et_id : et_uuids )
                {
                    absl::StrAppend(
                        &message, "; energy term uuid: ", et_id.string() );
                }
            }
            absl::StrAppend( &message, "\n\t - subset stat target list: " );

            for( const auto& [uuid, stats] : subset_stats_targets )
            {
                absl::StrAppend(
                    &message, "\n\t --> subset uuid: ", uuid.string() );
                for( const auto& stat : stats )
                {
                    absl::StrAppend( &message, "; stat: ", stat );
                }
            }
            absl::StrAppend(
                &message, "\n\t - subset sampler: ", subset_samplers.size() );
            for( const auto& [uuid, sampler] : subset_samplers )
            {
                absl::StrAppend(
                    &message, "\n\t --> subset uuid: ", uuid.string() );
            }
            absl::StrAppend( &message, "\n\t END " );
            return message;
        }
    };

    UserProblem create_problems(
        const MultitypeStraussDescription& description )
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
            auto subset_id = problem.object_set.add_subset();
            problem.object_subset_id.push_back( subset_id );

            problem.subset_samplers.emplace( subset_id,
                geode::UniformPointSetSampler< 2 >{ problem.box, subset_id } );
            OPENGEODE_EXCEPTION( points_desc.death_birth_ratio > 0.,
                "Object cannot be add or removed. Please set a BIRTH-DEATH "
                "with a positive probability." );
            proposal_kernel->add_move(
                std::make_unique< geode::BirthDeathMove< geode::Point2D > >(
                    problem.subset_samplers.at( subset_id ),
                    points_desc.death_birth_ratio, points_desc.birth_ratio ) );
            if( points_desc.change_ratio > 0. )
            {
                proposal_kernel->add_move(
                    std::make_unique< geode::ChangeMove< geode::Point2D > >(
                        problem.subset_samplers.at( subset_id ),
                        points_desc.change_ratio ) );
            }

            // energy terms here we define intra subset terms (can be
            // several of them) need to add inter set energy terms
            std::vector< geode::uuid > energy_terms;
            energy_terms.push_back( problem.gibbs_energy.add_energy_term(
                std::make_unique< geode::DensityTerm< geode::Point2D > >(
                    points_desc.name, points_desc.density, subset_id ) ) );

            auto interaction = std::make_unique<
                geode::EuclideanCutoffInteraction< geode::Point2D > >(
                points_desc.distance_treshold, points_desc.scope );
            energy_terms.push_back( problem.gibbs_energy.add_energy_term(
                std::make_unique< geode::PairwiseTerm< geode::Point2D > >(
                    "interaction", points_desc.gamma, std::move( interaction ),
                    subset_id ) ) );
            problem.subset_energy_term_ids.emplace( subset_id, energy_terms );

            std::vector< double > expected_statistics;
            expected_statistics.push_back(
                points_desc.expected_number_of_objects );
            expected_statistics.push_back(
                points_desc.expected_number_of_intersections );
            problem.subset_stats_targets.emplace(
                subset_id, expected_statistics );
        }
        DEBUG( proposal_kernel->string() );
        problem.mh_sampler =
            std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                problem.gibbs_energy, std::move( proposal_kernel ) );

        return problem;
    }
    void test_convergence(
        const MultitypeStraussDescription& problem_description,
        geode::RandomEngine& engine )
    {
        auto problem = create_problems( problem_description );
        SDEBUG( problem );

        problem.mh_sampler->walk( problem.object_set, engine, 500 );

        //        std::vector< double > sum_points(
        //        problem.object_subset_id.size(), 0. ); std::vector< double >
        //        sum_nb_interactions(
        //            problem.object_subset_id.size(), 0. );
        //
        //        auto N = problem_description.nb_realizations;

        //        for( const auto i : geode::Range{ N } )
        //        {
        //            problem.mh_sampler->walk(
        //                problem.object_set, engine,
        //                problem_description.nb_steps );
        //            for( const auto subset_id :
        //                geode::Range{ problem.object_subset_id.size() } )
        //            {
        //                const auto& energy_term_uuids =
        //                    problem.subset_energy_term_ids.at(
        //                        problem.object_subset_id[subset_id] );
        //                sum_points[subset_id] +=
        //                    problem.gibbs_energy.energy_term_statistic(
        //                        problem.object_set, energy_term_uuids[0] );
        //                sum_nb_interactions[subset_id] +=
        //                    problem.gibbs_energy.energy_term_statistic(
        //                        problem.object_set, energy_term_uuids[1] );
        //            }
        //        }
        //        std::transform( sum_points.begin(), sum_points.end(),
        //            sum_points.begin(), [N]( double p ) {
        //                return p / N;
        //            } );
        //        std::transform( sum_nb_interactions.begin(),
        //        sum_nb_interactions.end(),
        //            sum_nb_interactions.begin(), [N]( double p ) {
        //                return p / N;
        //            } );
        //        for( const auto subset_id :
        //            geode::Range{ problem.object_subset_id.size() } )
        //        {
        //            const auto& subset_uuid =
        //            problem.object_subset_id[subset_id]; const auto&
        //            expected_stats =
        //                problem.subset_stats_targets.at( subset_uuid );
        //
        //            geode::Logger::info( "[MH test] mean points = ",
        //                sum_points[subset_id], " (expected ",
        //                expected_stats[0],
        //                ")"
        //                " and mean interactions = ",
        //                sum_nb_interactions[subset_id], " (expected ",
        //                expected_stats[1], ")" );
        //            const auto error_stat =
        //                std::abs( sum_points[subset_id] - expected_stats[0] )
        //                / expected_stats[0];
        //            //            OPENGEODE_EXCEPTION( error_stat < 0.015,
        //            //                "[MH test] mean number of points not
        //            close to
        //            //                enought to " " expected value-- > error:
        //            ",
        //            //                error_stat );
        //            if( expected_stats[1] == 0 )
        //            {
        //                OPENGEODE_EXCEPTION(
        //                    sum_nb_interactions[subset_id] <
        //                    geode::GLOBAL_EPSILON,
        //                    "[MH test] Number of interactions not close to
        //                    enought to " " expected value-- > error : ",
        //                    sum_nb_interactions[subset_id] );
        //            }
        //            else
        //            {
        //                const auto error_interactions =
        //                    std::abs(
        //                        sum_nb_interactions[subset_id] -
        //                        expected_stats[1] )
        //                    / expected_stats[1];
        //                //                OPENGEODE_EXCEPTION(
        //                error_interactions < 0.1,
        //                //                    "[MH test] Number of
        //                interactions not
        //                //                    close to enought to " " expected
        //                value-- >
        //                //                    error : ", error_interactions );
        //            }
        //        }
        SDEBUG( problem );
    }

    void test_single_type_poisson()
    {
        geode::Logger::info( "TEST - MH SINGLE TYPE STRAUSS" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-Strauss-single@" );
        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        std::array< double, 5 > nb_points{ 22.6, 27.4, 31.3, 36.1, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4, 8, 13, 36 };

        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            MultitypeStraussDescription problem_description;
            problem_description.min_point = geode::Point2D{ { 0., 0. } };
            problem_description.max_point = geode::Point2D{ { 10., 10. } };

            StraussDescription description1;
            description1.name = "set1";
            description1.density = 0.5;
            description1.expected_number_of_objects = nb_points[config];

            description1.scope =
                geode::PairwiseInteraction< geode::Point2D >::SCOPE::all_subset;
            description1.distance_treshold = 1;
            description1.gamma = gamma_values[config];
            description1.expected_number_of_intersections =
                nb_interactions[config];

            // mh dynamic
            description1.death_birth_ratio = 1.;
            description1.birth_ratio = 0.5;
            description1.change_ratio = 0.;

            std::vector< StraussDescription > poisson_description{
                description1
            };
            problem_description.set_desc = poisson_description;

            problem_description.nb_steps = 1000.;
            problem_description.nb_realizations = 1000.;

            test_convergence( problem_description, engine );
        }
        geode::Logger::info( "TEST - MH SINGLE TYPE STRAUSS ... SUCCESS!" );
    }

    void test_multi_type_strauss()
    {
        geode::Logger::info( "TEST - MH MULTI TYPE STRAUSS" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-Strauss-multi@" );
        std::array< double, 4 > gamma_values{ 0, 0.3, 0.7, 1.0 };
        std::array< double, 4 > nb_points{ 22.6, 27.4, 31.3, 36.1 };
        std::array< double, 4 > nb_interactions{ 0, 4, 8, 13 };

        MultitypeStraussDescription problem_description;
        problem_description.min_point = geode::Point2D{ { 0., 0. } };
        problem_description.max_point = geode::Point2D{ { 10., 10. } };

        StraussDescription description1;
        description1.name = "set1";
        description1.density = 0.5;
        description1.expected_number_of_objects = nb_points[0];

        description1.scope =
            geode::PairwiseInteraction< geode::Point2D >::SCOPE::same_subset;
        description1.distance_treshold = 1;
        description1.gamma = gamma_values[0];
        description1.expected_number_of_intersections = nb_interactions[0];

        // mh dynamic
        description1.death_birth_ratio = 1.;
        description1.birth_ratio = 0.5;
        description1.change_ratio = 0.;

        StraussDescription description2;
        description2.name = "set2";
        description2.density = 0.5;
        description2.expected_number_of_objects = nb_points[1];

        description2.scope =
            geode::PairwiseInteraction< geode::Point2D >::SCOPE::same_subset;
        description2.distance_treshold = 2;
        description2.gamma = gamma_values[1];
        description2.expected_number_of_intersections = nb_interactions[1];

        // mh dynamic
        description2.death_birth_ratio = 1.;
        description2.birth_ratio = 0.5;
        description2.change_ratio = 0.;

        StraussDescription description3;
        description3.name = "set3";
        description3.density = 0.3;
        description3.expected_number_of_objects = nb_points[2];

        description3.scope =
            geode::PairwiseInteraction< geode::Point2D >::SCOPE::same_subset;
        description3.distance_treshold = 1;
        description3.gamma = gamma_values[2];
        description3.expected_number_of_intersections = nb_interactions[2];

        // mh dynamic
        description3.death_birth_ratio = 1.;
        description3.birth_ratio = 0.5;
        description3.change_ratio = 0.;

        StraussDescription description4;
        description4.name = "set4";
        description4.density = 0.5;
        description4.expected_number_of_objects = nb_points[3];

        description4.scope =
            geode::PairwiseInteraction< geode::Point2D >::SCOPE::same_subset;
        description4.distance_treshold = 1;
        description4.gamma = gamma_values[3];
        description4.expected_number_of_intersections = nb_interactions[3];

        // mh dynamic
        description4.death_birth_ratio = 1.;
        description4.birth_ratio = 0.5;
        description4.change_ratio = 0.;

        std::vector< StraussDescription > strauss_description{ // description1,
            // description2,
            description3, description4
        };
        problem_description.set_desc = strauss_description;

        problem_description.nb_steps = 10000.;
        problem_description.nb_realizations = 1000.;

        test_convergence( problem_description, engine );

        geode::Logger::info( "TEST - MH SINGLE TYPE STRAUSS ... SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        // test_single_type_poisson();
        test_multi_type_strauss();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}