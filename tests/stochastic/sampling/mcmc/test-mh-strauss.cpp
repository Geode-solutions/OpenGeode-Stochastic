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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/models/energy_terms/density_term.hpp>
#include <geode/stochastic/models/energy_terms/gibbs_energy.hpp>
#include <geode/stochastic/models/energy_terms/pairwise_term.hpp>
#include <geode/stochastic/models/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/models/proposal/classical_proposals.hpp>
#include <geode/stochastic/models/simulation_runner.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

namespace
{
    struct SetDescription
    {
        std::string name;
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    struct PoissonDensityDescription
    {
        std::string name;
        double density;
        double target_count;
    };

    struct PairwiseInteractionDescription
    {
        std::vector< std::string > names;
        double strength;
        double distance_threshold;
        // geode::PairwiseInteraction::SCOPE interaction_scope;
        double target_interaction_count;
    };

    class StraussSimulationRunner
        : public geode::SimulationRunner< geode::Point2D >
    {
    public:
        StraussSimulationRunner( const geode::SpatialDomain< 2 >& domain )
            : geode::SimulationRunner< geode::Point2D >( domain )
        {
        }

        void add_set_descriptor( const SetDescription& descriptor )
        {
            set_descriptors_.push_back( descriptor );
        }

        void add_density_descriptor(
            const PoissonDensityDescription& descriptor )
        {
            density_descriptors_.push_back( descriptor );
        }

        void add_interaction_descriptor(
            const PairwiseInteractionDescription& descriptor )
        {
            interaction_descriptors_.push_back( descriptor );
        }

        void initialize() override
        {
            auto proposal_kernel =
                std::make_unique< geode::ProposalKernel< geode::Point2D > >();

            // Mapping set names -> UUID
            std::unordered_map< std::string, geode::uuid > name_to_uuid;

            // Step 1: create object sets and samplers
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = this->object_sets_.add_set( set_desc.name );
                name_to_uuid[set_desc.name] = set_id;

                this->set_samplers_.push_back(
                    std::make_unique< geode::UniformPointSetSampler< 2 > >(
                        this->domain_ ) );

                geode::add_birth_death_change_moves( this->set_samplers_.back(),
                    *proposal_kernel, set_id, set_desc.birth_ratio,
                    set_desc.death_ratio, set_desc.change_ratio );
            }

            // Step 2: create density energy terms
            for( const auto& density_desc : density_descriptors_ )
            {
                const auto set_id = name_to_uuid.at( density_desc.name );
                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique<
                            geode::DensityTerm< geode::Point2D > >(
                            absl::StrCat( density_desc.name, "_density" ),
                            density_desc.density,
                            std::vector< geode::uuid >{ set_id },
                            this->domain_ ) ) );

                this->ordered_target_statistics_.push_back(
                    density_desc.target_count );
            }

            // Step 3: create pairwise interaction terms
            for( const auto& interaction_desc : interaction_descriptors_ )
            {
                std::vector< geode::uuid > set_ids;
                for( const auto& name : interaction_desc.names )
                {
                    set_ids.emplace_back( name_to_uuid.at( name ) );
                }

                auto interaction = std::make_unique<
                    geode::EuclideanCutoffInteraction< geode::Point2D > >(
                    interaction_desc.distance_threshold
                    /*,interaction_desc.interaction_scope*/ );

                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique<
                            geode::PairwiseTerm< geode::Point2D > >(
                            absl::StrCat(
                                absl::StrJoin( interaction_desc.names, "_" ),
                                "_interaction" ),
                            interaction_desc.strength, set_ids,
                            std::move( interaction ), this->domain_ ) ) );

                this->ordered_target_statistics_.push_back(
                    interaction_desc.target_interaction_count );
            }

            this->mh_sampler_ =
                std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                    this->energy_terms_collection_,
                    std::move( proposal_kernel ) );
        }

        void check_statistics(
            const geode::StatisticsMonitor& statistic_monitoring ) const
        {
            const auto& computed_means = statistic_monitoring.means();

            for( const auto stat_id :
                geode::Range{ this->energy_terms_collection_.size() } )
            {
                const auto& term = energy_terms_collection_.get(
                    ordered_energy_terms_[stat_id] );

                const auto expected_mean =
                    this->ordered_target_statistics_[stat_id];
                auto target_vs_mean_error =
                    std::abs( computed_means[stat_id] - expected_mean );
                if( expected_mean > 0 )
                {
                    target_vs_mean_error /= expected_mean;
                }

                OPENGEODE_EXCEPTION( target_vs_mean_error < 0.1,
                    "[MH test] Statistic value ", computed_means[stat_id],
                    " for energy term: ",
                    term.name().value_or( term.id().string() ),
                    " not close enough to expected value ", expected_mean,
                    " --> error: ", target_vs_mean_error );
            }
        }

    private:
        std::vector< SetDescription > set_descriptors_;
        std::vector< PoissonDensityDescription > density_descriptors_;
        std::vector< PairwiseInteractionDescription > interaction_descriptors_;
    };

    void test_single_type_strauss()
    {
        geode::Logger::info(
            "TEST - MH SINGLE TYPE STRAUSS (with intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-STRAUSS@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 10.0, 10.0 } } );
        // todo change!!
        geode::SpatialDomain domain( box, 1. );

        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        std::array< double, 5 > nb_points{ 19.5, 24.4, 31.3, 36.1, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4, 8, 15, 43 };

        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            // --- Object set
            SetDescription setA;
            setA.name = "A";
            setA.birth_ratio = 1.0;
            setA.death_ratio = 1.0;
            setA.change_ratio = 1.0;

            // --- Density term
            PoissonDensityDescription densityA;
            densityA.name = "A";
            densityA.density = 0.5;
            densityA.target_count = nb_points[config];

            // --- Intra-set pairwise interaction (Strauss process)
            PairwiseInteractionDescription intraA;
            intraA.names = { "A" }; // same set
            intraA.strength = gamma_values[config];
            intraA.distance_threshold = 1;
            // intraA.interaction_scope =
            // geode::PairwiseInteraction::SCOPE::INTRA;
            intraA.target_interaction_count = nb_interactions[config];

            StraussSimulationRunner runner( domain );
            runner.add_set_descriptor( setA );
            runner.add_density_descriptor( densityA );
            runner.add_interaction_descriptor( intraA );

            runner.initialize();

            // run simulation
            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_strauss_test_", config );

            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 1000;
            sim_config.metropolis_hasting_steps = 1000;
            sim_config.burn_in_steps = 1000;
            sim_config.printer = printer_config;

            auto statistic_monitoring = runner.run( engine, sim_config );
            runner.check_statistics( statistic_monitoring );
        }

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_multitype_strauss()
    {
        geode::Logger::info(
            "TEST - MH MULTITYPE STRAUSS (with inter-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-multi-STRAUSS@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 10.0, 10.0 } } );
        // todo change!!
        geode::SpatialDomain domain( box, 0. );

        std::array< double, 3 > gamma_values{ 0, 0.5, 1.0 };
        std::array< double, 3 > nb_points01{ 3.5, 5, 10.0 };
        std::array< double, 3 > nb_points02{ 14, 21, 40.0 };
        std::array< double, 3 > nb_points03{ 11, 16, 30. };
        std::array< double, 3 > nb_interactions01{ 0, 15, 95 };
        std::array< double, 3 > nb_interactions02{ 8, 20, 85 };

        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            // --- Sets
            SetDescription set01{ "set01", 1.0, 3.0, 1.0 };
            SetDescription set02{ "set02", 3.0, 0.5, 1.0 };
            SetDescription set03{ "set03", 4.0, 1.0, 1.0 };

            // --- Density terms
            PoissonDensityDescription d01{ "set01", 0.1, nb_points01[config] };
            PoissonDensityDescription d02{ "set02", 0.4, nb_points02[config] };
            PoissonDensityDescription d03{ "set03", 0.3, nb_points03[config] };

            // --- Pairwise interactions
            // 1. Intra-type (repulsion within same set)
            PairwiseInteractionDescription intra01{ { "set01", "set02",
                                                        "set03" },
                gamma_values[config], 1., nb_interactions01[config] };
            PairwiseInteractionDescription intra02{ { "set02" }, 1., 2.,
                nb_interactions02[config] };

            StraussSimulationRunner runner( domain );
            runner.add_set_descriptor( set01 );
            runner.add_set_descriptor( set02 );
            runner.add_set_descriptor( set03 );

            runner.add_density_descriptor( d01 );
            runner.add_density_descriptor( d02 );
            runner.add_density_descriptor( d03 );

            runner.add_interaction_descriptor( intra01 );
            runner.add_interaction_descriptor( intra02 );

            runner.initialize();

            // run simulation
            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_multitype_strauss_test" );

            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 750;
            sim_config.metropolis_hasting_steps = 1000;
            sim_config.burn_in_steps = 1000;
            sim_config.printer = printer_config;

            auto statistic_monitoring = runner.run( engine, sim_config );
            runner.check_statistics( statistic_monitoring );
        }

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_single_type_strauss();
        test_multitype_strauss();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}