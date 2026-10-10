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
    constexpr double DOMAIN_SIZE{ 10. };

    // Poisson process (gamma = 1) of intensity lambda: expected number of
    // pairs closer than r whose middle is in the domain V (the buffer is
    // larger than r / 2): lambda^2 |V| pi r^2 / 2
    double expected_poisson_pairs( double lambda, double r )
    {
        return lambda * lambda * DOMAIN_SIZE * DOMAIN_SIZE * M_PI * r * r / 2.;
    }

    void test_single_type_strauss()
    {
        geode::Logger::info(
            "TEST - MH SINGLE TYPE STRAUSS (with intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-STRAUSS@" );

        // NOLINTBEGIN(*-magic-numbers)
        std::array< double, 5 > gamma_values{ 0, 0.3, 0.5, 0.7, 1.0 };
        // gamma < 1: reference values from long simulations
        // gamma = 1: Poisson process, exact values
        std::array< double, 5 > nb_points{ 21.5, 26.3, 30.2, 35.6, 50. };
        std::array< double, 5 > nb_interactions{ 0, 4.1, 8.1, 14.4,
            expected_poisson_pairs( 0.5, 1. ) };
        for( const auto config : geode::Range{ gamma_values.size() } )
        {
            geode::SpatialDomainConfig< 2 > domain;
            domain.min_point = geode::Point2D{ { 0, 0 } };
            domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
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
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        domain.buffer_size = 2.;

        std::array< double, 3 > gamma_values{ 0, 0.5, 1.0 };
        // gamma < 1: reference values from long simulations
        // gamma = 1: Poisson processes, exact values
        std::array< double, 3 > nb_points01{ 7.7, 8.7, 10.0 };
        std::array< double, 3 > nb_points02{ 19.2, 26.0, 40.0 };
        std::array< double, 3 > nb_points03{ 16.5, 21.1, 30. };
        std::array< double, 3 > nb_interactions01{ 0, 10.4,
            expected_poisson_pairs( 0.1, 1. )
                + expected_poisson_pairs( 0.4, 1. )
                + expected_poisson_pairs( 0.3, 1. ) };
        // set02 pairs (gamma = 1) are only modified by the first interaction
        std::array< double, 3 > nb_interactions02{ 18.8, 38.4,
            expected_poisson_pairs( 0.4, 2. ) };
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

    void test_attractive_strauss_is_rejected()
    {
        geode::Logger::info( "TEST - STRAUSS PROCESS WITH GAMMA > 1" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::StraussProcessBuilder< geode::Point2D > strauss;
        auto& set_config = strauss.add_set( "set_A", 0.5 );
        geode_unused( set_config );
        bool rejected{ false };
        try
        {
            strauss.add_interaction( { "set_A" }, 2., 1.0, std::nullopt );
        }
        catch( const geode::OpenGeodeException& )
        {
            rejected = true;
        }
        // NOLINTEND(*-magic-numbers)
        geode::OpenGeodeStochasticStochasticException::test( rejected,
            "[Strauss] attractive interaction (gamma > 1) should be "
            "rejected" );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_same_sets_interactions_have_unique_names()
    {
        geode::Logger::info(
            "TEST - STRAUSS INTERACTIONS ON SAME SETS HAVE UNIQUE NAMES" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0, 0 } };
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        domain.buffer_size = 2.;

        geode::StraussProcessBuilder< geode::Point2D > strauss;
        strauss.set_domain( domain );
        auto& set_config = strauss.add_set( "set_A", 0.5 );
        geode_unused( set_config );
        strauss.add_interaction( { "set_A" }, 0.5, 1.0, 3. );
        strauss.add_interaction( { "set_A" }, 0.5, 2.0, 10. );
        // NOLINTEND(*-magic-numbers)

        const auto& stats = strauss.expected_statistics();
        geode::OpenGeodeStochasticStochasticException::test(
            stats.size() == 2 && stats[0].term_name == "pwint_set_A"
                && stats[1].term_name == "pwint_set_A_2",
            "[Strauss] interactions on the same sets should get distinct "
            "names" );

        const auto context = strauss.build_simulation_context();
        geode::TargetStatistics target_stats{ *context.model, stats };
        geode::OpenGeodeStochasticStochasticException::test(
            target_stats.active_terms().size() == 2,
            "[Strauss] both interaction targets should be active" );

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
        test_attractive_strauss_is_rejected();
        test_same_sets_interactions_have_unique_names();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}