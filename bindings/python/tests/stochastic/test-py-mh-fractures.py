# -*- coding: utf-8 -*-
# Copyright (c) 2019 - 2025 Geode-solutions
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

    box = og.BoundingBox2D()
    box.add_point(og.Point2D([0.0, 0.0]))
    box.add_point(og.Point2D([100.0, 100.0]))
    domain = stochastic.SpatialDomain2D(box,10.)

    # --- Object set
    setA = stochastic.FractureSetDescription()
    setA.name = "A"
    setA.add_observed_fracture(og.Point2D([10.0, 10.0]), og.Point2D([20.0, 20.0]))

    # length
    setA.length.distribution_type =stochastic.DistributionType("UniformClosed")
    setA.length.min_value = 1.0
    setA.length.max_value = 10.0

    # azimuth
    setA.azimuth.distribution_type =stochastic.DistributionType("UniformClosed")
    setA.azimuth.min_value = 1.0
    setA.azimuth.max_value = 10.0

    # positioning
    setA.p20 = 0.06
    setA.minimal_spacing = 1.0

    runner = stochastic.FractureSimulationRunner(domain)
    runner.add_fracture_set_descriptor(setA)
    runner.initialize()

    # Simulation printer
    printer_config = stochastic.SimulationPrinterConfigurator()
    printer_config.output_folder = os.path.join(printer_config.output_folder , "py_single_fracture_set")

    sim_config = stochastic.SimulationConfigurator()
    sim_config.realizations = 500
    sim_config.metropolis_hasting_steps = 1000
    sim_config.burn_in_steps = 1000
    sim_config.printer = printer_config

    # Run simulation
    statistic_monitoring = runner.run(engine, sim_config)
    runner.check_statistics(statistic_monitoring)

    print("--> SUCCESS!")


def test_two_fracture_sets_simulator():
    print("TEST - MH TWO SET FRACTURE SIMULATOR (with intra-set interactions)")

    engine = stochastic.RandomEngine()
    engine.set_seed("@mh-test-single-Fracture-set@")

    box = og.BoundingBox2D()
    box.add_point(og.Point2D([0.0, 0.0]))
    box.add_point(og.Point2D([100.0, 100.0]))
    domain = stochastic.SpatialDomain2D(box,10.)
    # --- Object set A
    setA = stochastic.FractureSetDescription()
    setA.name = "A"
    setA.length.distribution_type = stochastic.DistributionType("UniformClosed")
    setA.length.min_value = 1.0
    setA.length.max_value = 10.0
    setA.azimuth.distribution_type = stochastic.DistributionType("VonMises")
    setA.azimuth.mean = 45
    setA.azimuth.kappa = 1.0
    setA.p20 = 0.05
    setA.minimal_spacing = 1.0

    # --- Object set B
    setB = stochastic.FractureSetDescription()
    setB.name = "B"
    setB.length.distribution_type =stochastic.DistributionType("TruncatedLogNormal")
    setB.length.min_value = 1.0
    setB.length.max_value = 50.0
    setB.length.mean = 1.0
    setB.length.standard_deviation = 1.0
    setB.azimuth.distribution_type =stochastic.DistributionType("UniformClosed")
    setB.azimuth.min_value = 90.0
    setB.azimuth.max_value = 100.0
    setB.p20 = 0.03
    setB.minimal_spacing = 2.0

    runner = stochastic.FractureSimulationRunner(domain)
    runner.add_fracture_set_descriptor(setA)
    runner.add_fracture_set_descriptor(setB)
    runner.initialize()

    printer_config = stochastic.SimulationPrinterConfigurator()
    printer_config.output_folder = os.path.join(printer_config.output_folder ,printer_config.output_folder, "py_two_fracture_sets")

    sim_config = stochastic.SimulationConfigurator()
    sim_config.realizations = 500
    sim_config.metropolis_hasting_steps = 1000
    sim_config.burn_in_steps = 1000
    sim_config.printer = printer_config

    statistic_monitoring = runner.run(engine, sim_config)
    runner.check_statistics(statistic_monitoring)

    print("--> SUCCESS!")


if __name__ == "__main__":

    test_fracture_simulator()
    test_two_fracture_sets_simulator()
