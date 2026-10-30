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

import geode_test_env

import opengeode as og
import opengeode_stochastic_py_stochastic as stochastic

DOMAIN_SIZE = 20.0
MIN_LENGTH = 0.5
MAX_LENGTH = 1.0


def fracture_domain():
    domain = stochastic.SpatialDomainConfig2D()
    domain.min_point = og.Point2D([0.0, 0.0])
    domain.max_point = og.Point2D([DOMAIN_SIZE, DOMAIN_SIZE])
    # buffer larger than the longest fracture: no edge effect
    domain.buffer_size = 2.0 * MAX_LENGTH
    return domain


def set_uniform_sampler(fset):
    fset.sampler.length.distribution_type = stochastic.DistributionType("UniformClosed")
    fset.sampler.length.min_value = MIN_LENGTH
    fset.sampler.length.max_value = MAX_LENGTH

    fset.sampler.azimuth.distribution_type = stochastic.DistributionType("UniformClosed")
    fset.sampler.azimuth.min_value = 0.0
    fset.sampler.azimuth.max_value = 180.0


def run_and_validate(engine, fractures, output_name):
    runner = fractures.build_simulation_runner()

    sim_config = stochastic.SimulationConfigurator()
    sim_config.realizations = 1000
    sim_config.metropolis_hasting_steps = 100
    sim_config.burn_in_steps = 1000

    printer_config = stochastic.SimulationPrinterConfigurator()
    printer_config.output_folder = os.path.join(printer_config.output_folder, output_name)
    sim_config.printer = printer_config

    statistic_tracker = runner.run(engine, sim_config)
    runner.validate_statistics(statistic_tracker, fractures.expected_statistics())
    return runner, statistic_tracker


def check_named_results(runner, tracker, set_names, term_names):
    if runner.term_names() != term_names or tracker.term_names() != term_names:
        raise ValueError("[Test] Wrong statistic names: " + str(tracker.term_names()))
    means = tracker.means_by_term()
    if means != dict(zip(term_names, tracker.means())):
        raise ValueError("[Test] Wrong means by term")
    if sorted(tracker.variances_by_term()) != sorted(term_names):
        raise ValueError("[Test] Wrong variances by term")

    state = runner.state_realization()
    if state.set_names() != set_names:
        raise ValueError("[Test] Wrong set names: " + str(state.set_names()))
    nb_segments = 0
    for set_name in set_names:
        segments = state.segments(set_name)
        if len(segments) != state.nb_objects_in_set(set_name) or any(len(s) != 4 for s in segments):
            raise ValueError("[Test] Wrong segments list size")
        lengths = [((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5 for x0, y0, x1, y1 in segments]
        nb_fixed = state.nb_fixed_objects_in_set(set_name)
        if any(length < MIN_LENGTH - 1e-6 or length > MAX_LENGTH + 1e-6 for length in lengths[nb_fixed:]):
            raise ValueError("[Test] Sampled segment length out of bounds")
        nb_segments += len(segments)
    if nb_segments != state.nb_objects():
        raise ValueError("[Test] Missing segments")


def test_fracture_set():
    print("TEST - FRACTURE SET P20 + P21")

    engine = stochastic.RandomEngine()
    engine.set_seed("@py-mh-test-fracture-set@")

    p20 = 0.05
    expected_count = p20 * DOMAIN_SIZE * DOMAIN_SIZE
    expected_length = expected_count * (MIN_LENGTH + MAX_LENGTH) / 2.0

    fractures = stochastic.FractureProcessBuilder()
    fractures.set_domain(fracture_domain())
    fset = fractures.add_fracture_set("fset_A", p20, expected_count)
    set_uniform_sampler(fset)
    # p21 = 1: no effect on the model, used to monitor the total length
    fractures.add_intensity("fset_A", 1.0, expected_length)

    runner, tracker = run_and_validate(engine, fractures, "py_fracture_set")
    check_named_results(runner, tracker, ["fset_A"], ["fset_A_p20", "fset_A_p21"])
    print("--> SUCCESS!")


def test_two_fracture_sets():
    print("TEST - TWO FRACTURE SETS (spacing, observations and X-nodes)")

    engine = stochastic.RandomEngine()
    engine.set_seed("@py-mh-test-two-fracture-sets@")

    fractures = stochastic.FractureProcessBuilder()
    fractures.set_domain(fracture_domain())

    fset_a = fractures.add_fracture_set("fset_A", 0.05)
    fset_b = fractures.add_fracture_set("fset_B", 0.03)
    # definitions stay valid after adding other sets
    set_uniform_sampler(fset_a)
    set_uniform_sampler(fset_b)
    fset_b.dynamics.change_ratio = 2.0
    if fset_b.dynamics.change_ratio != 2.0:
        raise ValueError("[Test] Set dynamics not modified in place")
    fset_b.sampler.azimuth.distribution_type = stochastic.DistributionType("VonMises")
    fset_b.sampler.azimuth.mean = 60.0
    fset_b.sampler.azimuth.kappa = 1.0

    fset_a.add_fixed_segment(og.Point2D([1.0, 11.0]), og.Point2D([11.0, 20.0]))
    if fset_a.nb_fixed_objects() != 1:
        raise ValueError("[Test] Wrong number of observed fractures")

    # hard-core: no fractures closer than the spacing
    fractures.add_minimal_spacing("fset_A", 0.5, 0.0)
    # inhibited X-nodes: sets never intersect each other
    fractures.add_x_node_interaction(["fset_A", "fset_B"], 0.0, 0.0)

    runner, tracker = run_and_validate(engine, fractures, "py_two_fracture_sets")
    check_named_results(
        runner,
        tracker,
        ["fset_A", "fset_B"],
        ["fset_A_p20", "fset_B_p20", "fset_A_spacing", "x_node_fset_A_fset_B"],
    )
    if runner.state_realization().nb_fixed_objects_in_set("fset_A") != 1:
        raise ValueError("[Test] Observed fracture missing in the realization")
    print("--> SUCCESS!")


def test_keyword_constructors():
    print("TEST - KEYWORD CONSTRUCTORS")

    domain = stochastic.SpatialDomainConfig2D(
        min_point=og.Point2D([0.0, 0.0]),
        max_point=og.Point2D([DOMAIN_SIZE, DOMAIN_SIZE]),
        buffer_size=2.0 * MAX_LENGTH,
    )
    if domain.buffer_size != 2.0 * MAX_LENGTH or domain.max_point.value(0) != DOMAIN_SIZE:
        raise ValueError("[Test] Wrong domain configuration")

    length = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=MIN_LENGTH, max_value=MAX_LENGTH
    )
    if length.min_value != MIN_LENGTH or length.mean is not None:
        raise ValueError("[Test] Wrong distribution description")
    if not length.distribution_type.matches(stochastic.DistributionType("UniformClosed")):
        raise ValueError("[Test] Wrong distribution type")

    fractures = stochastic.FractureProcessBuilder()
    fractures.set_domain(domain)
    fset = fractures.add_fracture_set("fset_A", 0.05)
    fset.sampler.length = length
    fset.sampler.azimuth = stochastic.DoubleDistributionDescription("VonMises", mean=30.0, kappa=10.0)
    # the sampler distributions keep their name when replaced
    if fset.sampler.length.name != "length" or fset.sampler.azimuth.name != "azimuth":
        raise ValueError("[Test] Sampler distribution names lost")
    if fset.sampler.azimuth.kappa != 10.0:
        raise ValueError("[Test] Sampler distribution not modified in place")
    # distribution types can be given as strings
    fset.sampler.azimuth.distribution_type = "UniformClosed"
    fset.sampler.azimuth.min_value = 0.0
    fset.sampler.azimuth.max_value = 180.0
    fset.dynamics = stochastic.ObjectSetDynamicsConfig(change_ratio=2.0)
    if fset.dynamics.change_ratio != 2.0 or fset.dynamics.birth_ratio != 1.0:
        raise ValueError("[Test] Wrong dynamics configuration")

    sampler = stochastic.Segment2DSamplerConfig(length=length, move_ratio=0.2)
    if sampler.length.name != "length" or sampler.move_ratio != 0.2 or sampler.rotation_ratio != 1.0:
        raise ValueError("[Test] Wrong sampler configuration")

    target = stochastic.TargetStatisticConfig("fset_A_p20", 20.0, tolerance=0.2)
    if (target.term_name, target.value, target.tolerance) != ("fset_A_p20", 20.0, 0.2):
        raise ValueError("[Test] Wrong target statistic")

    printer = stochastic.SimulationPrinterConfigurator(
        output_folder=os.path.join(os.getcwd(), "py_keyword_constructors"),
        realisations_print_frequency=50,
    )
    if printer.statistics_filename != "statistics.txt":
        raise ValueError("[Test] Printer defaults lost")
    sim_config = stochastic.SimulationConfigurator(
        realizations=100, metropolis_hasting_steps=100, burn_in_steps=1000, printer=printer
    )
    if sim_config.printer.realisations_print_frequency != 50:
        raise ValueError("[Test] Printer not set")
    if stochastic.SimulationConfigurator().printer is not None:
        raise ValueError("[Test] Default configurator should not print")

    engine = stochastic.RandomEngine()
    engine.set_seed("@py-keyword-constructors@")
    runner = fractures.build_simulation_runner()
    tracker = runner.run(engine, sim_config)
    runner.validate_statistics(tracker, [target])
    print("--> SUCCESS!")


if __name__ == "__main__":
    test_fracture_set()
    test_two_fracture_sets()
    test_keyword_constructors()
