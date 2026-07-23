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
    void test_single_type_strauss()
    {
        geode::Logger::info(
            "TEST - MH SINGLE TYPE STRAUSS (with intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-STRAUSS@" );

        // NOLINTBEGIN(*-magic-numbers)
        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        std::array< double, 5 > nb_points{ 20.8, 25.8, 29.9, 35.5, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4.7, 9.2, 16.5, 42.8 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            geode::SpatialDomainConfig< 2 > domain;
            domain.min_point = geode::Point2D{ { 0, 0 } };
            domain.max_point = geode::Point2D{ { 10, 10 } };
            domain.buffer_size = 2.;

            geode::StraussProcessBuilder< geode::Point2D > strauss;

            strauss.set_domain( domain );

            auto& set_config =
                strauss.add_set( "set_A", 0.5, nb_points[config] );
            strauss.add_interaction( { "set_A" }, gamma_values[config], 1.0,
                nb_interactions[config], true );

            geode::SimulationRunner< geode::Point2D > runner{
                strauss.build_simulation_context()
            };

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
            geode::TargetStatistics target_stats{ runner.model(),
                strauss.expected_statistics() };
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
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0, 0 } };
        domain.max_point = geode::Point2D{ { 10, 10 } };
        domain.buffer_size = 2.;

        std::array< double, 3 > gamma_values{ 0, 0.5, 1.0 };
        std::array< double, 3 > nb_points01{ 7.5, 8.6, 10.0 };
        std::array< double, 3 > nb_points02{ 18.5, 25.4, 40.0 };
        std::array< double, 3 > nb_points03{ 16.0, 21.0, 30. };
        std::array< double, 3 > nb_interactions01{ 0, 11.5, 43.2 };
        std::array< double, 3 > nb_interactions02{ 26.3, 49.5, 116.7 };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            geode::StraussProcessBuilder< geode::Point2D > strauss;
            strauss.set_domain( domain );
            auto& set_config =
                strauss.add_set( "set01", 0.1, nb_points01[config] );
            geode_unused( set_config );

            auto& set_config_02 =
                strauss.add_set( "set02", 0.4, nb_points02[config] );
            geode_unused( set_config_02 );

            auto& set_config_03 =
                strauss.add_set( "set03", 0.3, nb_points03[config] );
            geode_unused( set_config_03 );

            strauss.add_interaction( { "set01", "set02", "set03" },
                gamma_values[config], 1.0, nb_interactions01[config], true );
            strauss.add_interaction(
                { "set02" }, 1.0, 2.0, nb_interactions02[config], true );

            // --- Pairwise interactions
            // 1. Intra-type (repulsion within same set)

            // run simulation
            geode::SimulationRunner< geode::Point2D > runner{
                strauss.build_simulation_context()
            };

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
            geode::TargetStatistics target_stats{ runner.model(),
                strauss.expected_statistics() };
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
        test_multitype_strauss();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}