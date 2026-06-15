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

#include <geode/stochastic/applications/fractures.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>

namespace
{
    void test_fracture_simulator()
    {
        geode::Logger::info( "TEST - MH SINGLE SET FRACTURE SIMULATOR( with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-Fracture-set@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::FractureNetworkDescription fnet_desc;
        fnet_desc.fnet_name = "One_Set_FNet";

        fnet_desc.domain.min_point = geode::Point2D{ { 0, 0 } };
        fnet_desc.domain.max_point = geode::Point2D{ { 100.0, 100.0 } };
        fnet_desc.domain.buffer_size = 10.;

        // --- Object set
        auto& fset = fnet_desc.add_fracture_set( "fset_A" );
        fset.sampler.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.length.min_value = 1;
        fset.sampler.length.max_value = 10.;

        fset.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.azimuth.min_value = 1;
        fset.sampler.azimuth.max_value = 10.;

        fset.p20 = 0.05;
        fset.p21 = 200;
        fset.minimal_spacing = 1.;

        // observed fractures
        fset.observed_fractures.push_back( { geode::Point2D{ { 0.0, 15. } },
            geode::Point2D{ { 15., 15. } } } );
        fset.observed_fractures.push_back( { geode::Point2D{ { 1.0, 11. } },
            geode::Point2D{ { 11., 20. } } } );

        geode::Logger::info( fnet_desc.string() );

        // runner
        auto context = build_fractures_simulation_context( fnet_desc );
        geode::FractureSimulationRunner runner{ std::move( context ) };
        // run simulation
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 500;
        sim_config.metropolis_hasting_steps = 100;
        sim_config.burn_in_steps = 1000;

        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/sim_one_fracture_set_test" );
        sim_config.printer = printer_config;

        auto statistic_tracker = runner.run( engine, sim_config );

        // NOLINTEND(*-magic-numbers)

        //        const auto targeted_statistics_descriptors =
        //            build_fractures_targeted_stat( fnet_desc );
        //        geode::TargetStatistics target_stats{ runner.model(),
        //            targeted_statistics_descriptors };
        //        geode::statistics::validate( statistic_tracker, target_stats
        //        );
        //
        const auto& fset_state = runner.state_realization().get_set(
            runner.state_realization().get_set_uuid(
                fnet_desc.fracture_sets[0].fset_name ) );

        geode::OpenGeodeStochasticStochasticException::test(
            fset_state.nb_fixed_objects() == 2,
            "nd fixed object = ", fset_state.nb_fixed_objects() );
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_two_fracture_sets_simulator()
    {
        geode::Logger::info( "TEST - MH TWO SET FRACTURE SIMULATOR (with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-two-Fracture-set@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::FractureNetworkDescription fnet_desc;
        fnet_desc.fnet_name = "Two_Set_FNet";
        constexpr double DOMAIN_BUFFER{ 10 };
        fnet_desc.domain = { geode::Point2D{ { 0, 0 } },
            geode::Point2D{ { 100.0, 100.0 } }, DOMAIN_BUFFER };

        // --- Object set
        auto& fset_01 = fnet_desc.add_fracture_set( "fset_01" );
        fset_01.sampler.length.distribution_type =
            geode::TruncatedPowerLaw::distribution_type_static();
        fset_01.sampler.length.alpha = 2.;
        fset_01.sampler.length.min_value = 1;
        fset_01.sampler.length.max_value = 10.;

        fset_01.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset_01.sampler.azimuth.min_value = 1;
        fset_01.sampler.azimuth.max_value = 10.;

        fset_01.p20 = 0.05;
        fset_01.p21 = 200;
        fset_01.minimal_spacing = 1.;

        auto& fset_02 = fnet_desc.add_fracture_set( "fset_02" );
        fset_02.sampler.length.distribution_type =
            geode::TruncatedLogNormal::distribution_type_static();
        fset_02.sampler.length.mean = 1;
        fset_02.sampler.length.standard_deviation = 1.;
        fset_02.sampler.length.min_value = 1;
        fset_02.sampler.length.max_value = 50.;

        fset_02.sampler.azimuth.distribution_type =
            geode::VonMises::distribution_type_static();
        fset_02.sampler.azimuth.mean = 60.;
        fset_02.sampler.azimuth.kappa = 1.;

        fset_02.p20 = 0.05;
        fset_02.p21 = 200;
        fset_02.minimal_spacing = 2.;

        fnet_desc.add_x_node_monitoring( 0.3 );

        // runner
        auto context = build_fractures_simulation_context( fnet_desc );
        geode::FractureSimulationRunner runner{ std::move( context ) };
        // run simulation
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 500;
        sim_config.metropolis_hasting_steps = 100;
        sim_config.burn_in_steps = 1000;

        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/sim_one_fracture_set_test" );
        sim_config.printer = printer_config;

        auto statistic_tracker = runner.run( engine, sim_config );

        // NOLINTEND(*-magic-numbers)

        //        const auto targeted_statistics_descriptors =
        //            build_fractures_targeted_stat( fnet_desc );
        //        geode::TargetStatistics target_stats{ runner.model(),
        //            targeted_statistics_descriptors };
        //        geode::statistics::validate( statistic_tracker, target_stats
        //        );
        //
        const auto& fset_state = runner.state_realization().get_set(
            runner.state_realization().get_set_uuid(
                fnet_desc.fracture_sets[0].fset_name ) );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_fracture_simulator();
        test_two_fracture_sets_simulator();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}