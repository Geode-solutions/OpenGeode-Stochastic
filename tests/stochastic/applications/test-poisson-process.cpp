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

#include <geode/stochastic/applications/poisson_process.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>

#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace
{

    void test_single_type_poisson()
    {
        geode::Logger::info( "TEST - MH SINGLE TYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-POISSON@" );

        // NOLINTBEGIN(*-magic-numbers)
        std::array< double, 4 > birth_ratio{ 0.1, 0.5, 2., 4. };
        std::array< double, 4 > change_ratio{ 0., 1., 1., 0. };

        for( const auto config : geode::Range{ birth_ratio.size() } )
        {
            geode::PoissonProcessDescription< geode::Point2D > poisson;
            poisson.domain = { geode::Point2D{ { 0, 0 } },
                geode::Point2D{ { 10, 10 } }, 0. };
            auto& set_config = poisson.add_set( "set_A", "density_A" );
            set_config.lambda = 0.3;
            set_config.expected_nb_objects = 30;
            set_config.birth_ratio = birth_ratio[config];
            set_config.death_ratio = 1.0;
            set_config.change_ratio = change_ratio[config];

            auto simulation_context = build_poisson_process( poisson );
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
                    "/sim_point_poisson_test_", config );
            sim_config.printer = printer_config;

            auto statistic_tracker = runner.run( engine, sim_config );

            const auto targeted_statistics_descriptors =
                build_poisson_targeted_stat( poisson );
            geode::TargetStatistics target_stats{ runner.model(),
                targeted_statistics_descriptors };

            geode::statistics::validate( statistic_tracker, target_stats );
        }
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_multitype_poisson()
    {
        geode::Logger::info( "TEST - MH MULTITYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-POISSON-multi@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::PoissonProcessDescription< geode::Point2D > poisson;
        poisson.domain = { geode::Point2D{ { 0, 0 } },
            geode::Point2D{ { 10, 10 } }, 0. };
        auto& set_config_01 = poisson.add_set( "set01", "density_01" );
        set_config_01.lambda = 0.1;
        set_config_01.expected_nb_objects = 10;
        set_config_01.birth_ratio = 2.0;
        set_config_01.death_ratio = 3.0;
        set_config_01.change_ratio = 1.0;

        auto& set_config_02 = poisson.add_set( "set02", "density_02" );
        set_config_02.lambda = 0.4;
        set_config_01.expected_nb_objects = 40;
        set_config_02.birth_ratio = 3.0;
        set_config_02.death_ratio = 0.5;
        set_config_02.change_ratio = 1.0;

        auto& set_config_03 = poisson.add_set( "set03", "density_03" );
        set_config_03.lambda = 0.3;
        set_config_01.expected_nb_objects = 30;
        set_config_03.birth_ratio = 4.0;
        set_config_03.death_ratio = 1.0;
        set_config_03.change_ratio = 1.0;

        auto context = build_poisson_process( poisson );
        geode::SimulationRunner< geode::Point2D > runner{ std::move(
            context ) };

        // run simulation
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 2000;
        sim_config.metropolis_hasting_steps = 100;
        sim_config.burn_in_steps = 1000;

        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/sim_point_multitype_poisson_test" );
        sim_config.printer = printer_config;

        auto statistic_tracker = runner.run( engine, sim_config );

        // NOLINTEND(*-magic-numbers)

        const auto targeted_statistics_descriptors =
            build_poisson_targeted_stat( poisson );
        geode::TargetStatistics target_stats{ runner.model(),
            targeted_statistics_descriptors };
        geode::statistics::validate( statistic_tracker, target_stats );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_single_type_poisson();
        test_multitype_poisson();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}