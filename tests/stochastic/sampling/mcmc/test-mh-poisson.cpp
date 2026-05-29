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

#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>
#include <geode/stochastic/models/gibbs_energy.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

#include <geode/stochastic/sampling/mcmc/proposal/object_set_dynamic_config.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>

namespace
{

    using PoissonDensityDescription = geode::SingleObjectTermConfig;

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
            geode::SimulationContextConfig< geode::Point2D > poisson_config;

            poisson_config.domain = { geode::Point2D{ { 0.0, 0.0 } },
                geode::Point2D{ { 10.0, 10.0 } }, 0. };
            // --- Set description
            geode::ObjectSetDefinition< geode::Point2D > set_a;
            set_a.name = "A";
            set_a.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
            set_a.dynamics.birth_ratio = birth_ratio[config];
            set_a.dynamics.death_ratio = 1.0;
            set_a.dynamics.change_ratio = change_ratio[config];
            poisson_config.sets.emplace_back( set_a );

            // --- Energy term description
            PoissonDensityDescription density_a;
            density_a.term_name = "density";
            density_a.object_set_names = { "A" };
            density_a.lambda = 0.3;
            density_a.object_feature = geode::ObjectInDomainFeatureConfig{};
            poisson_config.model.terms.emplace_back( density_a );

            auto context = build_simulation_context( poisson_config );
            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                context ) };

            // run simulation

            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 1000;
            sim_config.metropolis_hasting_steps = 1000;
            sim_config.burn_in_steps = 1000;

            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_poisson_test_", config );
            sim_config.printer = printer_config;

            auto statistic_tracker = runner.run( engine, sim_config );

            std::vector< geode::TargetStatisticConfig >
                targeted_statistics_descriptors;
            geode::TargetStatisticConfig stat_a{ "density", 30.0, 0.15 };
            targeted_statistics_descriptors.push_back( stat_a );

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

        geode::SimulationContextConfig< geode::Point2D > poisson_config;

        // NOLINTBEGIN(*-magic-numbers)
        poisson_config.domain = { geode::Point2D{ { 0.0, 0.0 } },
            geode::Point2D{ { 10.0, 10.0 } }, 0. };

        // --- Set descriptions
        geode::ObjectSetDefinition< geode::Point2D > set01;
        set01.name = "set01";
        set01.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
        set01.dynamics.birth_ratio = 2.0;
        set01.dynamics.death_ratio = 3.0;
        set01.dynamics.change_ratio = 1.0;
        poisson_config.sets.emplace_back( set01 );

        geode::ObjectSetDefinition< geode::Point2D > set02;
        set02.name = "set02";
        set02.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
        set02.dynamics.birth_ratio = 3.0;
        set02.dynamics.death_ratio = 0.5;
        set02.dynamics.change_ratio = 1.0;
        poisson_config.sets.emplace_back( set02 );

        geode::ObjectSetDefinition< geode::Point2D > set03;
        set03.name = "set03";
        set03.sampler = geode::ObjectSamplerConfig< geode::Point2D >{};
        set03.dynamics.birth_ratio = 4.0;
        set03.dynamics.death_ratio = 1.0;
        set03.dynamics.change_ratio = 1.0;
        poisson_config.sets.emplace_back( set03 );

        // --- Energy term descriptions
        PoissonDensityDescription density01;
        density01.term_name = "density01";
        density01.object_set_names = { "set01" };
        density01.lambda = 0.1;
        density01.object_feature = geode::ObjectInDomainFeatureConfig{};
        poisson_config.model.terms.emplace_back( density01 );

        PoissonDensityDescription density02;
        density02.term_name = "density02";
        density02.object_set_names = { "set02" };
        density02.lambda = 0.4;
        density02.object_feature = geode::ObjectInDomainFeatureConfig{};
        poisson_config.model.terms.emplace_back( density02 );

        PoissonDensityDescription density03;
        density03.term_name = "density03";
        density03.object_set_names = { "set03" };
        density03.lambda = 0.3;
        density03.object_feature = geode::ObjectInDomainFeatureConfig{};
        poisson_config.model.terms.emplace_back( density03 );

        auto context = build_simulation_context( poisson_config );
        geode::SimulationRunner< geode::Point2D > runner{ std::move(
            context ) };

        // run simulation
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 1500;
        sim_config.metropolis_hasting_steps = 1000;
        sim_config.burn_in_steps = 3000;

        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/sim_point_multitype_poisson_test" );
        sim_config.printer = printer_config;

        auto statistic_tracker = runner.run( engine, sim_config );

        std::vector< geode::TargetStatisticConfig >
            targeted_statistics_descriptors;
        geode::TargetStatisticConfig stat01{ "density01", 10.0, 0.15 };
        targeted_statistics_descriptors.push_back( stat01 );
        geode::TargetStatisticConfig stat02{ "density02", 40.0, 0.15 };
        targeted_statistics_descriptors.push_back( stat02 );
        geode::TargetStatisticConfig stat03{ "density03", 30.0, 0.15 };
        targeted_statistics_descriptors.push_back( stat03 );
        // NOLINTEND(*-magic-numbers)

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