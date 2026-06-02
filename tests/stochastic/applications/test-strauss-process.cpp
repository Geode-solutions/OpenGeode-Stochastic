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

#include <geode/stochastic/applications/strauss_process.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>

#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

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
            geode::StraussProcessDescription< geode::Point2D > strauss;

            strauss.domain = { geode::Point2D{ { 0.0, 0.0 } },
                geode::Point2D{ { 10.0, 10.0 } }, 1. };

            auto& set_config = strauss.add_set( "set_A", "density_A" );
            set_config.lambda = 0.5;
            set_config.expected_nb_objects = nb_points[config];

            auto& interaction_config =
                strauss.add_interaction( "interaction_A" );
            interaction_config.set_names = { "set_A" };
            interaction_config.gamma = gamma_values[config];
            interaction_config.distance = 1.0;
            interaction_config.expected_nb_interactions =
                nb_interactions[config];

            auto simulation_context = geode::build_strauss_process( strauss );
            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                simulation_context ) };

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

            const auto targeted_statistics_descriptors =
                geode::build_strauss_targeted_stat( strauss );

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
            geode::StraussProcessDescription< geode::Point2D > strauss;
            strauss.domain = { geode::Point2D{ { 0, 0 } },
                geode::Point2D{ { 10, 10 } }, 2. };

            auto& set_config_01 = strauss.add_set( "set01", "density_01" );
            set_config_01.lambda = 0.1;
            set_config_01.expected_nb_objects = 10;
            set_config_01.birth_ratio = 1.0;
            set_config_01.death_ratio = 3.0;
            set_config_01.change_ratio = 1.0;

            auto& set_config_02 = strauss.add_set( "set02", "density_02" );
            set_config_02.lambda = 0.4;
            set_config_01.expected_nb_objects = 40;
            set_config_02.birth_ratio = 3.0;
            set_config_02.death_ratio = 0.5;
            set_config_02.change_ratio = 1.0;

            auto& set_config_03 = strauss.add_set( "set03", "density_03" );
            set_config_03.lambda = 0.3;
            set_config_01.expected_nb_objects = 30;
            set_config_03.birth_ratio = 4.0;
            set_config_03.death_ratio = 1.0;
            set_config_03.change_ratio = 1.0;

            auto& interaction_config_01 =
                strauss.add_interaction( "interaction_01" );
            interaction_config_01.set_names = { "set01", "set02", "set03" };
            interaction_config_01.gamma = gamma_values[config];
            interaction_config_01.distance = 1.0;
            interaction_config_01.expected_nb_interactions =
                nb_interactions01[config];

            auto& interaction_config_02 =
                strauss.add_interaction( "interaction_02" );
            interaction_config_02.set_names = { "set02" };
            interaction_config_02.gamma = 1.;
            interaction_config_02.distance = 2.0;
            interaction_config_02.expected_nb_interactions =
                nb_interactions02[config];

            // --- Pairwise interactions
            // 1. Intra-type (repulsion within same set)

            // run simulation
            auto context = build_strauss_process( strauss );
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

            const auto targeted_statistics_descriptors =
                geode::build_strauss_targeted_stat( strauss );
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
        test_single_type_strauss();
        // test_multitype_strauss();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}