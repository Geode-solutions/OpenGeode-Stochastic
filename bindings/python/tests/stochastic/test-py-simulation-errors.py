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

import gc

import geode_test_env

import opengeode as og
import opengeode_stochastic_py_stochastic as stochastic


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def raises(error, function, message):
    try:
        function()
    except error as exception:
        return exception
    raise ValueError("[Test] " + message)


def domain_config(size=10.0, buffer_size=2.0):
    return stochastic.SpatialDomainConfig2D(
        min_point=og.Point2D([0.0, 0.0]),
        max_point=og.Point2D([size, size]),
        buffer_size=buffer_size,
    )


def small_builder(p20=0.05):
    builder = stochastic.FractureProcessBuilder()
    builder.set_domain(domain_config())
    fset = builder.add_fracture_set("fset", p20, 5.0)
    fset.sampler.length = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.5, max_value=1.0
    )
    fset.sampler.azimuth = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.0, max_value=3.14
    )
    return builder, fset


def test_exception_hierarchy():
    # the exception of this repository reaches Python as the OpenGeode hierarchy
    engine = stochastic.RandomEngine()
    exception = raises(
        og.OpenGeodeDataError,
        lambda: engine.sample_bernoulli(2.0),
        "No data error",
    )
    check(isinstance(exception, og.OpenGeodeError), "Not an OpenGeodeError")
    check(isinstance(exception, RuntimeError), "Not a RuntimeError")
    check(exception.type == "data", "Wrong type: %r" % exception.type)
    check("Bernoulli" in exception.message, "Wrong message: %r" % exception.message)
    check(exception.project == "OpenGeodeStochastic", "Wrong project: %r" % exception.project)
    check(exception.library, "No library")
    check(exception.message in exception.details, "Details do not hold the message")
    print("test_exception_hierarchy: OK")


def test_domain_errors():
    box = og.BoundingBox2D()
    box.add_point(og.Point2D([0.0, 0.0]))
    raises(og.OpenGeodeDataError, lambda: stochastic.SpatialDomain2D(box, 1.0), "Empty box")
    box.add_point(og.Point2D([2.0, 4.0]))
    raises(og.OpenGeodeDataError, lambda: stochastic.SpatialDomain2D(box, -1.0), "Negative buffer")
    domain = stochastic.SpatialDomain2D(box, 0.5)
    check(domain.buffer_size() == 0.5, "buffer_size")
    check(domain.n_volume() == 8.0, "n_volume")
    check(domain.extended_n_volume() == 3.0 * 5.0, "extended_n_volume")
    check(domain.contains(og.Point2D([1.0, 1.0])), "contains")
    check(not domain.contains(og.Point2D([-0.25, 1.0])), "contains outside")
    check(domain.extended_contains(og.Point2D([-0.25, 1.0])), "extended_contains")
    check(not domain.extended_contains(og.Point2D([-1.0, 1.0])), "extended_contains outside")
    check(domain.box().max().value(1) == 4.0, "box")
    check(domain.extended_box().max().value(1) == 4.5, "extended_box")
    check("buffer" in domain.string(), "string")
    box3 = og.BoundingBox3D()
    box3.add_point(og.Point3D([0.0, 0.0, 0.0]))
    box3.add_point(og.Point3D([1.0, 1.0, 1.0]))
    check(stochastic.SpatialDomain3D(box3, 0.0).n_volume() == 1.0, "3D volume")
    # the references keep the domain alive
    reference = stochastic.SpatialDomain2D(box, 0.0).box()
    gc.collect()
    check(reference.max().value(0) == 2.0, "The box does not keep its domain alive")
    print("test_domain_errors: OK")


def test_builder_errors():
    builder = stochastic.FractureProcessBuilder()
    # no domain: the default configuration is empty
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "Empty domain accepted")

    builder, _ = small_builder()
    raises(
        og.OpenGeodeDataError,
        lambda: builder.add_x_node_interaction(["fset"], 0.5),
        "X-node with one set",
    )
    raises(
        og.OpenGeodeDataError,
        lambda: builder.add_x_node_interaction(["fset", "other"], 1.5),
        "beta > 1",
    )
    raises(
        og.OpenGeodeDataError,
        lambda: builder.add_x_node_interaction(["fset", "other"], -0.5),
        "beta < 0",
    )
    # sets of the terms are looked up when the simulation is built
    builder.add_intensity("unknown_set", 1.0)
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "Unknown set accepted")

    # an unusable sampler: no length distribution type
    builder = stochastic.FractureProcessBuilder()
    builder.set_domain(domain_config())
    builder.add_fracture_set("fset", 0.1)
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "Sampler without laws")

    # a length larger than the buffer
    builder, fset = small_builder()
    fset.sampler.length.max_value = 5.0
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "Length above the buffer")

    # no change move at all
    builder, fset = small_builder()
    fset.sampler.translation_ratio = 0.0
    fset.sampler.rotation_ratio = 0.0
    fset.sampler.stretch_ratio = 0.0
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "No change move")
    builder, fset = small_builder()
    fset.sampler.rotation_ratio = -1.0
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "Negative change move")
    builder, fset = small_builder()
    fset.dynamics = stochastic.ObjectSetDynamicsConfig(birth_ratio=0.0, death_ratio=0.0, change_ratio=0.0)
    raises(og.OpenGeodeDataError, lambda: builder.build_simulation_runner(), "No move at all")
    print("test_builder_errors: OK")


def test_object_sets_errors():
    runner = small_builder()[0].build_simulation_runner()
    sets = runner.state_realization()
    check(sets.set_names() == ["fset"], "set_names")
    check(sets.nb_sets() == 1 and sets.nb_objects() == 0, "Empty state")
    for function in (sets.nb_objects_in_set, sets.nb_fixed_objects_in_set, sets.segments):
        raises(og.OpenGeodeDataError, lambda: function("unknown"), "Unknown set accepted")
    check(sets.segments("fset") == [], "No segment before the run")
    print("test_object_sets_errors: OK")


def test_validation():
    engine = stochastic.RandomEngine()
    engine.set_seed(3)
    builder, _ = small_builder()
    runner = builder.build_simulation_runner()
    config = stochastic.SimulationConfigurator(realizations=10, metropolis_hasting_steps=20, burn_in_steps=50)
    tracker = runner.run(engine, config)
    check(tracker.statiscal_count() == 10, "Wrong count: %d" % tracker.statiscal_count())
    # the expected value is far from the mean: result error
    far = [stochastic.TargetStatisticConfig("fset_p20", 1.0e6)]
    exception = raises(og.OpenGeodeResultError, lambda: runner.validate_statistics(tracker, far), "No result error")
    check(exception.type == "result", "Wrong type: " + exception.type)
    check("fset_p20" in exception.message, "The message does not name the term")
    # a very large tolerance validates
    runner.validate_statistics(tracker, [stochastic.TargetStatisticConfig("fset_p20", 1.0e6, tolerance=10.0)])
    # the expected statistics of the builder are valid targets
    check(
        [target.term_name for target in builder.expected_statistics()] == ["fset_p20"],
        "Wrong expected statistics",
    )
    # unknown term
    raises(
        og.OpenGeodeDataError,
        lambda: runner.validate_statistics(tracker, [stochastic.TargetStatisticConfig("nope", 1.0)]),
        "Unknown term accepted",
    )
    # a tracker of another runner has other terms
    other = builder.build_simulation_runner()
    other_tracker = other.run(engine, config)
    raises(
        og.OpenGeodeDataError,
        lambda: runner.validate_statistics(other_tracker, [stochastic.TargetStatisticConfig("fset_p20", 5.0, tolerance=100.0)]),
        "A tracker of another runner was accepted",
    )
    print("test_validation: OK")


if __name__ == "__main__":
    test_exception_hierarchy()
    test_domain_errors()
    test_builder_errors()
    test_object_sets_errors()
    test_validation()
