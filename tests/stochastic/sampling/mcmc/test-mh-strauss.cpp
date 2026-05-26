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
#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>
#include <geode/stochastic/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
// #include <geode/stochastic/spatial/pairwise_interactions.hpp>
#include <geode/stochastic/inference/statistics_tools.hpp>

namespace
{
    struct SetDescription
    {
        std::string name;
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    using PoissonDensityDescription = geode::SingleObjectTermConfig;
    using PairwiseInteractionDescription = geode::PairwiseTermConfig;

    class StraussSimulationRunner
        : public geode::SimulationRunner< geode::Point2D >
    {
    public:
        explicit StraussSimulationRunner(
            const geode::SpatialDomain< 2 >& domain )
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

        void add_target_statistics(
            const geode::TargetStatisticConfig& statistic_descriptor )
        {
            targeted_statistics_descriptors_.push_back( statistic_descriptor );
        }
        std::unique_ptr< geode::ProposalKernel< geode::Point2D > >
            create_sets_and_set_samplers()
        {
            auto proposal_kernel =
                std::make_unique< geode::ProposalKernel< geode::Point2D > >();

            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = this->object_sets_.add_set( set_desc.name );
                this->set_samplers_.push_back(
                    std::make_unique< geode::UniformPointSetSampler< 2 > >(
                        domain_ ) );

                geode::add_birth_death_change_moves( this->set_samplers_.back(),
                    *proposal_kernel, set_id, set_desc.birth_ratio,
                    set_desc.death_ratio, set_desc.change_ratio );
            }
            return proposal_kernel;
        }

        void create_model()
        {
            geode::ModelConfig config;
            for( const auto& energy_desc : density_descriptors_ )
            {
                config.terms.emplace_back( energy_desc );
            }
            for( const auto& interaction_desc : interaction_descriptors_ )
            {
                config.terms.emplace_back( interaction_desc );
            }

            model_ = std::move( geode::build_model< geode::Point2D >(
                config, object_sets_, domain_ ) );
            create_target_statistics();
        }

        void create_target_statistics()
        {
            target_statistics_.emplace( *model_ );
            for( const auto& target_stat : targeted_statistics_descriptors_ )
            {
                target_statistics_->set_target( target_stat );
            }
        }

        void initialize() override
        {
            auto proposal_kernel = create_sets_and_set_samplers();
            create_model();

            this->mh_sampler_ =
                std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                    *model_, std::move( proposal_kernel ) );
            create_target_statistics();
        }

    private:
        std::vector< SetDescription > set_descriptors_;
        std::vector< PoissonDensityDescription > density_descriptors_;
        std::vector< PairwiseInteractionDescription > interaction_descriptors_;

        std::vector< geode::TargetStatisticConfig >
            targeted_statistics_descriptors_;
    };

    void test_single_type_strauss()
    {
        geode::Logger::info(
            "TEST - MH SINGLE TYPE STRAUSS (with intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-STRAUSS@" );
        // NOLINTBEGIN(*-magic-numbers)
        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 10.0, 10.0 } } );
        geode::SpatialDomain domain( box, 1. );

        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        std::array< double, 5 > nb_points{ 19.5, 24.4, 31.3, 36.1, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4.7, 9.8, 18.7, 50.3 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            // --- Object set
            SetDescription set_a;
            set_a.name = "A";
            set_a.birth_ratio = 1.0;
            set_a.death_ratio = 1.0;
            set_a.change_ratio = 1.0;

            // --- Energy term description
            PoissonDensityDescription density_a;
            density_a.term_name = "density_a";
            density_a.object_set_names = { "A" };
            density_a.lambda = 0.5;
            density_a.object_feature = geode::ObjectInDomainFeatureConfig{};

            geode::TargetStatisticConfig stat_a{ "density_a", nb_points[config],
                0.1 };

            // --- Intra-set pairwise interaction (Strauss process)
            PairwiseInteractionDescription interaction_a;
            interaction_a.term_name = "interactionA";
            interaction_a.object_set_names_interactions = { { "A", "A" } };
            interaction_a.gamma = gamma_values[config];

            interaction_a.interaction_config =
                geode::MinimalDistanceCutoffConfig{
                    1.
                }; /// ici cela devrait etre un paramètre utilisateur

            geode::TargetStatisticConfig stat_intra_a{ "interactionA",
                nb_interactions[config], 0.1 };

            StraussSimulationRunner runner( domain );
            runner.add_set_descriptor( set_a );
            runner.add_density_descriptor( density_a );
            runner.add_interaction_descriptor( interaction_a );
            runner.add_target_statistics( stat_a );
            runner.add_target_statistics( stat_intra_a );
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

            auto statistic_tracker = runner.run( engine, sim_config );
            geode::statistics::validate(
                statistic_tracker, runner.target_statistics() );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_multitype_strauss()
    {
        geode::Logger::info(
            "TEST - MH MULTITYPE STRAUSS (with inter-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-multi-STRAUSS@" );
        // NOLINTBEGIN(*-magic-numbers)
        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 10.0, 10.0 } } );
        geode::SpatialDomain domain( box, 2. );

        std::array< double, 3 > gamma_values{ 0, 0.5, 1.0 };
        std::array< double, 3 > nb_points01{ 6.7, 8, 10.0 };
        std::array< double, 3 > nb_points02{ 17.5, 24.6, 40.0 };
        std::array< double, 3 > nb_points03{ 14.6, 19.4, 30. };
        std::array< double, 3 > nb_interactions01{ 0, 15, 59.8 };
        std::array< double, 3 > nb_interactions02{ 37.2, 70, 174 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            // --- Sets
            SetDescription set01{ "set01", 1.0, 3.0, 1.0 };
            SetDescription set02{ "set02", 3.0, 0.5, 1.0 };
            SetDescription set03{ "set03", 4.0, 1.0, 1.0 };

            // --- Density terms
            PoissonDensityDescription d01;
            d01.term_name = "density_set01";
            d01.object_set_names = { "set01" };
            d01.lambda = 0.1;
            d01.object_feature = geode::ObjectInDomainFeatureConfig{};

            geode::TargetStatisticConfig stat01{ "density_set01",
                nb_points01[config], 0.1 };

            PoissonDensityDescription d02;
            d02.term_name = "density_set02";
            d02.object_set_names = { "set02" };
            d02.lambda = 0.4;
            d02.object_feature = geode::ObjectInDomainFeatureConfig{};

            geode::TargetStatisticConfig stat02{ "density_set02",
                nb_points02[config], 0.1 };

            PoissonDensityDescription d03;
            d03.term_name = "density_set03";
            d03.object_set_names = { "set03" };
            d03.lambda = 0.3;
            d03.object_feature = geode::ObjectInDomainFeatureConfig{};

            geode::TargetStatisticConfig stat03{ "density_set03",
                nb_points03[config], 0.1 };

            // --- Pairwise interactions
            // 1. Intra-type (repulsion within same set)
            PairwiseInteractionDescription intra01;
            intra01.term_name = "interaction01";
            intra01.object_set_names_interactions = { { "set01", "set01" },
                { "set02", "set02" }, { "set03", "set03" } };
            intra01.gamma = gamma_values[config];

            intra01.interaction_config = geode::MinimalDistanceCutoffConfig{
                1.
            }; /// ici cela devrait etre un paramètre utilisateur

            geode::TargetStatisticConfig stat_intra_01{ "interaction01",
                nb_interactions01[config], 0.1 };

            PairwiseInteractionDescription intra02;
            intra02.term_name = "interaction02";
            intra02.object_set_names_interactions = { { "set02", "set02" } };
            intra02.gamma = 1.; // gamma_values[config];

            intra02.interaction_config = geode::MinimalDistanceCutoffConfig{
                2.
            }; /// ici cela devrait etre un paramètre utilisateur

            geode::TargetStatisticConfig stat_intra_02{ "interaction02",
                nb_interactions02[config], 0.1 };

            StraussSimulationRunner runner( domain );
            runner.add_set_descriptor( set01 );
            runner.add_set_descriptor( set02 );
            runner.add_set_descriptor( set03 );

            runner.add_density_descriptor( d01 );
            runner.add_density_descriptor( d02 );
            runner.add_density_descriptor( d03 );

            runner.add_interaction_descriptor( intra01 );
            runner.add_interaction_descriptor( intra02 );

            runner.add_target_statistics( stat01 );
            runner.add_target_statistics( stat02 );
            runner.add_target_statistics( stat03 );
            runner.add_target_statistics( stat_intra_01 );
            runner.add_target_statistics( stat_intra_02 );

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

            auto statistic_tracker = runner.run( engine, sim_config );
            geode::statistics::validate(
                statistic_tracker, runner.target_statistics() );
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