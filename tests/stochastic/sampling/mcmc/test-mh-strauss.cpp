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

#include <geode/stochastic/inference/statistics_tools.hpp>
#include <geode/stochastic/inference/target_statistics.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>

namespace
{
    using PoissonDensityDescription = geode::SingleObjectTermConfig;
    using PairwiseInteractionDescription = geode::PairwiseTermConfig;

    void test_single_type_strauss()
    {
        geode::Logger::info(
            "TEST - MH SINGLE TYPE STRAUSS (with intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-STRAUSS@" );

        // NOLINTBEGIN(*-magic-numbers)
        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        std::array< double, 5 > nb_points{ 19.5, 24.4, 31.3, 36.1, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4.7, 9.8, 18.7, 50.3 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            geode::SimulationContextConfig< geode::Point2D >
                simulation_cxt_config;

            simulation_cxt_config.domain = { geode::Point2D{ { 0.0, 0.0 } },
                geode::Point2D{ { 10.0, 10.0 } }, 1. };

            // --- Object set
            geode::ObjectSetDefinition< geode::Point2D > set_a;
            set_a.name = "A";
            set_a.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
            set_a.dynamics.birth_ratio = 1.0;
            set_a.dynamics.death_ratio = 1.0;
            set_a.dynamics.change_ratio = 1.0;
            simulation_cxt_config.sets.emplace_back( set_a );

            // --- Model - Energy term description
            PoissonDensityDescription density_a;
            density_a.term_name = "density_a";
            density_a.object_set_names = { "A" };
            density_a.lambda = 0.5;
            density_a.object_feature = geode::ObjectInDomainFeatureConfig{};
            simulation_cxt_config.model.terms.emplace_back( density_a );

            // --- Intra-set pairwise interaction (Strauss process)
            PairwiseInteractionDescription interaction_a;
            interaction_a.term_name = "interactionA";
            interaction_a.object_set_names_interactions = { { "A", "A" } };
            interaction_a.gamma = gamma_values[config];
            interaction_a.interaction_config =
                geode::MinimalDistanceCutoffConfig{
                    1.
                }; /// ici cela devrait etre un paramètre utilisateur
            simulation_cxt_config.model.terms.emplace_back( interaction_a );

            auto context = build_simulation_context( simulation_cxt_config );
            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                context ) };

            // run simulation
            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 2000;
            sim_config.metropolis_hasting_steps = 100;
            sim_config.burn_in_steps = 1000;

            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_strauss_test_", config );
            sim_config.printer = printer_config;

            auto statistic_tracker = runner.run( engine, sim_config );

            std::vector< geode::TargetStatisticConfig >
                targeted_statistics_descriptors;
            geode::TargetStatisticConfig stat_a{ "density_a", nb_points[config],
                0.1 };
            geode::TargetStatisticConfig stat_intra_a{ "interactionA",
                nb_interactions[config], 0.1 };
            targeted_statistics_descriptors.push_back( stat_a );
            targeted_statistics_descriptors.push_back( stat_intra_a );
            geode::TargetStatistics target_stats{ runner.model(),
                targeted_statistics_descriptors };

            geode::statistics::validate( statistic_tracker, target_stats );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_multitype_strauss()
    {
        geode::Logger::info(
            "TEST - MH MULTITYPE STRAUSS (with inter-set interactions ) " );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-multi-STRAUSS@" );

        // NOLINTBEGIN(*-magic-numbers)
        std::array< double, 3 > gamma_values{ 0, 0.5, 1.0 };
        std::array< double, 3 > nb_points01{ 6.7, 8, 10.0 };
        std::array< double, 3 > nb_points02{ 17.5, 24.6, 40.0 };
        std::array< double, 3 > nb_points03{ 14.6, 19.4, 30. };
        std::array< double, 3 > nb_interactions01{ 0, 15, 59.8 };
        std::array< double, 3 > nb_interactions02{ 37.2, 70, 174 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            geode::SimulationContextConfig< geode::Point2D >
                simulation_cxt_config;

            simulation_cxt_config.domain = { geode::Point2D{ { 0.0, 0.0 } },
                geode::Point2D{ { 10.0, 10.0 } }, 2. };

            // --- Object set
            geode::ObjectSetDefinition< geode::Point2D > set01;
            set01.name = "set01";
            set01.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
            set01.dynamics.birth_ratio = 1.0;
            set01.dynamics.death_ratio = 3.0;
            set01.dynamics.change_ratio = 1.0;
            simulation_cxt_config.sets.emplace_back( set01 );

            geode::ObjectSetDefinition< geode::Point2D > set02;
            set02.name = "set02";
            set02.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
            set02.dynamics.birth_ratio = 3.0;
            set02.dynamics.death_ratio = 0.5;
            set02.dynamics.change_ratio = 1.0;
            simulation_cxt_config.sets.emplace_back( set02 );

            geode::ObjectSetDefinition< geode::Point2D > set03;
            set03.name = "set03";
            set03.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
            set03.dynamics.birth_ratio = 4.0;
            set03.dynamics.death_ratio = 1.0;
            set03.dynamics.change_ratio = 1.0;
            simulation_cxt_config.sets.emplace_back( set03 );

            // --- Model - Energy term description
            // density
            PoissonDensityDescription density01;
            density01.term_name = "density01";
            density01.object_set_names = { "set01" };
            density01.lambda = 0.1;
            density01.object_feature = geode::ObjectInDomainFeatureConfig{};
            simulation_cxt_config.model.terms.emplace_back( density01 );

            PoissonDensityDescription density02;
            density02.term_name = "density02";
            density02.object_set_names = { "set02" };
            density02.lambda = 0.4;
            density02.object_feature = geode::ObjectInDomainFeatureConfig{};
            simulation_cxt_config.model.terms.emplace_back( density02 );

            PoissonDensityDescription density03;
            density03.term_name = "density03";
            density03.object_set_names = { "set03" };
            density03.lambda = 0.3;
            density03.object_feature = geode::ObjectInDomainFeatureConfig{};
            simulation_cxt_config.model.terms.emplace_back( density03 );

            // --- Pairwise interactions
            // 1. Intra-type (repulsion within same set)
            PairwiseInteractionDescription interaction_01;
            interaction_01.term_name = "interaction_01";
            interaction_01.object_set_names_interactions = {
                { "set01", "set01" }, { "set02", "set02" }, { "set03", "set03" }
            };
            interaction_01.gamma = gamma_values[config];
            interaction_01.interaction_config =
                geode::MinimalDistanceCutoffConfig{
                    1.
                }; /// ici cela devrait etre un paramètre utilisateur
            simulation_cxt_config.model.terms.emplace_back( interaction_01 );

            PairwiseInteractionDescription interaction_02;
            interaction_02.term_name = "interaction_02";
            interaction_02.object_set_names_interactions = { { "set02",
                "set02" } };
            interaction_02.gamma = 1.;
            interaction_02.interaction_config =
                geode::MinimalDistanceCutoffConfig{
                    2.
                }; /// ici cela devrait etre un paramètre utilisateur
            simulation_cxt_config.model.terms.emplace_back( interaction_02 );

            // run simulation
            auto context = build_simulation_context( simulation_cxt_config );
            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                context ) };

            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 2000;
            sim_config.metropolis_hasting_steps = 100;
            sim_config.burn_in_steps = 1000;

            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_multitype_strauss_test" );
            sim_config.printer = printer_config;

            auto statistic_tracker = runner.run( engine, sim_config );

            std::vector< geode::TargetStatisticConfig >
                targeted_statistics_descriptors;
            geode::TargetStatisticConfig stat01{ "density01",
                nb_points01[config], 0.1 };
            geode::TargetStatisticConfig stat02{ "density02",
                nb_points02[config], 0.1 };
            geode::TargetStatisticConfig stat03{ "density03",
                nb_points03[config], 0.1 };
            geode::TargetStatisticConfig stat_intra_01{ "interaction_01",
                nb_interactions01[config], 0.1 };
            geode::TargetStatisticConfig stat_intra_02{ "interaction_02",
                nb_interactions02[config], 0.1 };
            targeted_statistics_descriptors.push_back( stat01 );
            targeted_statistics_descriptors.push_back( stat02 );
            targeted_statistics_descriptors.push_back( stat03 );
            targeted_statistics_descriptors.push_back( stat_intra_01 );
            targeted_statistics_descriptors.push_back( stat_intra_02 );

            geode::TargetStatistics target_stats{ runner.model(),
                targeted_statistics_descriptors };

            geode::statistics::validate( statistic_tracker, target_stats );
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
        //  test_single_type_strauss();
        test_multitype_strauss();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}