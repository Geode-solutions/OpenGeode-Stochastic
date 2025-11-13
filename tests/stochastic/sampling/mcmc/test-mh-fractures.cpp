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

#include <geode/stochastic/sampling/mcmc/helpers/fracture_simulation_runner.hpp>

namespace
{
    void test_fracture_simulator()
    {
        geode::Logger::info( "TEST - MH SINGLE SET FRACTURE SIMULATOR (with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-Fracture-set@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 100.0, 100.0 } } );

        // --- Object set
        geode::FractureSetDescription setA;
        setA.name = "A";

        // length
        setA.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.length.min_value = 1;
        setA.length.max_value = 10.;

        // azimuth
        setA.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.azimuth.min_value = 1;
        setA.azimuth.max_value = 10.;

        // positionning
        setA.p20 = 0.05;
        setA.minimal_spacing = 1.;

        geode::FractureSimulationRunner runner( box );
        runner.add_fracture_set_descriptor( setA );

        runner.initialize();

        // run simulation
        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/single_fracture_set" );

        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 1000;
        sim_config.metropolis_hasting_steps = 1000;
        sim_config.burn_in_steps = 1000;
        sim_config.printer = printer_config;

        auto statistic_monitoring = runner.run( engine, sim_config );
        runner.check_statistics( statistic_monitoring );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_two_fracture_sets_simulator()
    {
        geode::Logger::info( "TEST - MH TWO SET FRACTURE SIMULATOR (with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-Fracture-set@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 100.0, 100.0 } } );

        // --- Object set
        geode::FractureSetDescription setA;
        setA.name = "A";

        // length
        setA.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.length.min_value = 1;
        setA.length.max_value = 10.;

        // azimuth
        setA.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.azimuth.min_value = 1;
        setA.azimuth.max_value = 10.;

        // positionning
        setA.p20 = 0.05;
        setA.minimal_spacing = 1.;

        // --- Object set
        geode::FractureSetDescription setB;
        setB.name = "B";

        // length
        setB.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setB.length.min_value = 1;
        setB.length.max_value = 10.;

        // azimuth
        setB.azimuth.distribution_type =
            geode::VonMises::distribution_type_static();
        setB.azimuth.mean = 60.;
        setB.azimuth.standard_deviation = 15.;

        // positionning
        setB.p20 = 0.05;
        setB.minimal_spacing = 2.;

        geode::FractureSimulationRunner runner( box );
        runner.add_fracture_set_descriptor( setA );
        runner.add_fracture_set_descriptor( setB );

        runner.initialize();

        // run simulation
        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder =
            absl::StrCat( printer_config.output_folder, "/two_fracture_sets" );

        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 1000;
        sim_config.metropolis_hasting_steps = 1000;
        sim_config.burn_in_steps = 1000;
        sim_config.printer = printer_config;

        auto statistic_monitoring = runner.run( engine, sim_config );
        runner.check_statistics( statistic_monitoring );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
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