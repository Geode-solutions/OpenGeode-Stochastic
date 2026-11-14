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
import threading
import time

import geode_test_env

import opengeode as og
import opengeode_stochastic_py_stochastic as stochastic


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def make_builder():
    builder = stochastic.FractureProcessBuilder()
    builder.set_domain(
        stochastic.SpatialDomainConfig2D(
            min_point=og.Point2D([0.0, 0.0]),
            max_point=og.Point2D([20.0, 20.0]),
            buffer_size=2.0,
        )
    )
    fset = builder.add_fracture_set("fset", 0.2)
    fset.sampler.length = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.5, max_value=1.0
    )
    fset.sampler.azimuth = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.0, max_value=3.14
    )
    return builder


def count_while_running(task):
    """Run task in a thread and count the loops of the main thread meanwhile"""
    errors = []

    def target():
        try:
            task()
        except BaseException as error:  # noqa: B036
            errors.append(error)

    thread = threading.Thread(target=target)
    ticks = 0
    thread.start()
    started = time.perf_counter()
    while thread.is_alive():
        ticks += 1
        time.sleep(0.001)
    elapsed = time.perf_counter() - started
    thread.join()
    if errors:
        raise errors[0]
    return ticks, elapsed


def test_run_releases_the_gil():
    engine = stochastic.RandomEngine()
    engine.set_seed(5)
    runner = make_builder().build_simulation_runner()
    steps = 200000

    ticks, elapsed = count_while_running(lambda: runner.run(engine, steps))
    # with the GIL held by the run, the main thread could not tick: it needs
    # about one tick per millisecond of the run
    check(elapsed > 0.05, "The run is too short to be measured: %f s" % elapsed)
    check(ticks > elapsed * 100, "The GIL is not released by run(engine, steps): %d ticks in %f s" % (ticks, elapsed))

    config = stochastic.SimulationConfigurator(
        realizations=20, metropolis_hasting_steps=steps // 20, burn_in_steps=steps
    )
    ticks, elapsed = count_while_running(lambda: runner.run(engine, config))
    check(elapsed > 0.05, "The simulation is too short to be measured: %f s" % elapsed)
    check(ticks > elapsed * 100, "The GIL is not released by run(engine, config): %d ticks in %f s" % (ticks, elapsed))
    print("test_run_releases_the_gil: OK")


def test_build_releases_the_gil():
    builder = make_builder()
    # the runners are built one after the other by two threads, the result must
    # be the same as in one thread
    runners = []
    threads = [threading.Thread(target=lambda: runners.append(builder.build_simulation_runner())) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    check(len(runners) == 4, "Runners are missing")
    check(all(runner.term_names() == ["fset_p20"] for runner in runners), "Wrong terms")
    print("test_build_releases_the_gil: OK")


def test_results_keep_the_runner_alive():
    engine = stochastic.RandomEngine()
    engine.set_seed(5)
    config = stochastic.SimulationConfigurator(realizations=5, metropolis_hasting_steps=10, burn_in_steps=10)

    runner = make_builder().build_simulation_runner()
    tracker = runner.run(engine, config)
    del runner
    gc.collect()
    check(tracker.term_names() == ["fset_p20"], "The tracker lost its runner")
    check(list(tracker.means_by_term()) == ["fset_p20"], "means_by_term")
    check(tracker.statiscal_count() == 5, "count")
    check(len(tracker.variances()) == 1 and len(tracker.means()) == 1, "variances and means")

    runner = make_builder().build_simulation_runner()
    sets = runner.run(engine, 50)
    del runner
    gc.collect()
    check(sets.nb_sets() == 1 and sets.nb_objects() == len(sets.segments("fset")), "The state lost its runner")

    runner = make_builder().build_simulation_runner()
    state = runner.state_realization()
    del runner
    gc.collect()
    check(state.set_names() == ["fset"], "state_realization lost its runner")

    builder = make_builder()
    definition = builder.add_fracture_set("second", 0.1)
    del builder
    gc.collect()
    definition.dynamics.birth_ratio = 3.0
    check(definition.name == "second" and definition.dynamics.birth_ratio == 3.0, "The set definition lost its builder")
    print("test_results_keep_the_runner_alive: OK")


def test_runner_is_independent_of_the_builder():
    engine = stochastic.RandomEngine()
    engine.set_seed(5)
    builder = make_builder()
    runner = builder.build_simulation_runner()
    builder.add_fracture_set("later", 0.1)
    check(runner.term_names() == ["fset_p20"], "A later change of the builder reached the runner")
    check(runner.state_realization().set_names() == ["fset"], "A later set reached the runner")
    print("test_runner_is_independent_of_the_builder: OK")


if __name__ == "__main__":
    test_run_releases_the_gil()
    test_build_releases_the_gil()
    test_results_keep_the_runner_alive()
    test_runner_is_independent_of_the_builder()
