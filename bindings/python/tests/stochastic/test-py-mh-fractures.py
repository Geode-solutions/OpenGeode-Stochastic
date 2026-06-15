# -*- coding: utf-8 -*-
# Copyright (c) 2019 - 2026 Geode-solutions
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import sys
import platform

if sys.version_info >= (3, 8, 0) and platform.system() == "Windows":
    for path in [x.strip() for x in os.environ["PATH"].split(";") if x]:
        os.add_dll_directory(path)

import opengeode as og
import opengeode_stochastic_py_stochastic as stochastic

from pathlib import Path

def test_fracture_simulator():
    print("TEST - MH SINGLE SET FRACTURE SIMULATOR (with intra-set interactions)")

    engine = stochastic.RandomEngine()
    engine.set_seed("@mh-test-single-Fracture-set@")

    fnet_desc = stochastic.FractureNetworkDescription()
    fnet_desc.name = "One_Set_FNet"

    fnet_desc.domain.min_point = og.Point2D([0.0, 0.0])
    fnet_desc.domain.max_point = og.Point2D([100.0, 100.0])
    fnet_desc.domain.buffer_size = 10.

    fset = fnet_desc.add_fracture_set( "fset_A" )
    
    fset.sampler.length.distribution_type = stochastic.DistributionType("UniformClosed")
    fset.sampler.length.min_value = 1
    fset.sampler.length.max_value = 10.

    fset.sampler.azimuth.distribution_type = stochastic.DistributionType("UniformClosed")
    fset.sampler.azimuth.min_value = 1.
    fset.sampler.azimuth.max_value = 10.

    fset.p20 = 0.05
    fset.p21 = 200
    fset.minimal_spacing = 1.

    fset.add_observed_fracture(og.Point2D([10.0, 10.0]), og.Point2D([20.0, 20.0]))
    fset.add_observed_fracture(og.Point2D([15.0, 15.0]), og.Point2D([0.0, 15.0]))
    fset.add_observed_fracture(og.Point2D([1.0, 11.0]), og.Point2D([11.0, 20.0]))

    print( fnet_desc.string() );

    runner = stochastic.build_fractures_simulation_runner( fnet_desc ) 

    sim_config = stochastic.SimulationConfigurator()
    sim_config.realizations = 500
    sim_config.metropolis_hasting_steps = 100
    sim_config.burn_in_steps = 1000

    printer_config = stochastic.SimulationPrinterConfigurator()
    printer_config.output_folder = os.path.join(printer_config.output_folder , "py_single_fracture_set")
    sim_config.printer = printer_config

    statistic_tracker = runner.run( engine, sim_config );
#    runner.check_statistics(statistic_monitoring)
    print("--> SUCCESS!")


def test_two_fracture_sets_simulator():
    print("TEST - MH TWO SET FRACTURE SIMULATOR (with intra-set interactions)")

    engine = stochastic.RandomEngine()
    engine.set_seed("@mh-test-two-Fracture-set@")

    fnet_desc = stochastic.FractureNetworkDescription()
    fnet_desc.name = "Two_Sets_FNet"

    fnet_desc.domain.min_point = og.Point2D([0.0, 0.0])
    fnet_desc.domain.max_point = og.Point2D([100.0, 100.0])
    fnet_desc.domain.buffer_size = 10.

    fset = fnet_desc.add_fracture_set( "fset_A" )
    
    fset.sampler.length.distribution_type = stochastic.DistributionType("UniformClosed")
    fset.sampler.length.min_value = 1
    fset.sampler.length.max_value = 10.

    fset.sampler.azimuth.distribution_type = stochastic.DistributionType("VonMises")
    fset.sampler.azimuth.mean = 45.
    fset.sampler.azimuth.kappa = 1.

    fset.p20 = 0.05
    fset.p21 = 200
    fset.minimal_spacing = 1.

    fset_b = fnet_desc.add_fracture_set( "fset_B" )
    
    fset_b.sampler.length.distribution_type = stochastic.DistributionType("TruncatedLogNormal")
    fset_b.sampler.length.min_value = 1
    fset_b.sampler.length.max_value = 50.
    fset_b.sampler.length.mean = 1.
    fset_b.sampler.length.standard_deviation = 1.0

    fset_b.sampler.azimuth.distribution_type =stochastic.DistributionType("UniformClosed")
    fset_b.sampler.azimuth.min_value = 90.0
    fset_b.sampler.azimuth.max_value = 100.0

    fset_b.p20 = 0.03
    fset_b.p21 = 100
    fset_b.minimal_spacing = 2.

    print( fnet_desc.string() );

    runner = stochastic.build_fractures_simulation_runner( fnet_desc ) 

    sim_config = stochastic.SimulationConfigurator()
    sim_config.realizations = 500
    sim_config.metropolis_hasting_steps = 100
    sim_config.burn_in_steps = 1000

    printer_config = stochastic.SimulationPrinterConfigurator()
    printer_config.output_folder = os.path.join(printer_config.output_folder , "py_two_fracture_set")
    sim_config.printer = printer_config

    statistic_tracker = runner.run( engine, sim_config );
#    runner.check_statistics(statistic_monitoring)
    print("--> SUCCESS!")

if __name__ == "__main__":

    test_fracture_simulator()
    test_two_fracture_sets_simulator()
