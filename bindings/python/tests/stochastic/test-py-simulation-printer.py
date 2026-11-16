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
import tempfile

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


def small_builder():
    builder = stochastic.FractureProcessBuilder()
    builder.set_domain(
        stochastic.SpatialDomainConfig2D(
            min_point=og.Point2D([0.0, 0.0]),
            max_point=og.Point2D([10.0, 10.0]),
            buffer_size=2.0,
        )
    )
    fset = builder.add_fracture_set("fset", 0.05)
    fset.sampler.length = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.5, max_value=1.0
    )
    fset.sampler.azimuth = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=0.0, max_value=3.14
    )
    return builder


def test_printer_errors():
    engine = stochastic.RandomEngine()
    engine.set_seed(3)
    builder = small_builder()
    runner = builder.build_simulation_runner()
    with tempfile.TemporaryDirectory() as folder:
        # a zero frequency used to divide by zero while printing
        printer = stochastic.SimulationPrinterConfigurator(
            output_folder=folder, realisations_print_frequency=0
        )
        config = stochastic.SimulationConfigurator(
            realizations=3, metropolis_hasting_steps=1, burn_in_steps=1, printer=printer
        )
        raises(og.OpenGeodeDataError, lambda: runner.run(engine, config), "Zero print frequency accepted")
        # with the realization printing off the frequency is not used
        printer.print_realisations = False
        # the configurator holds a copy of the printer
        check(config.printer.print_realisations, "The configurator shares the printer")
        config.printer = printer
        runner.run(engine, config)
        check(os.path.exists(os.path.join(folder, "statistics.txt")), "No statistics file")
        # an output folder that is a file
        blocker = os.path.join(folder, "blocker")
        open(blocker, "w").close()
        printer = stochastic.SimulationPrinterConfigurator(
            output_folder=os.path.join(blocker, "sub"), print_realisations=False
        )
        config.printer = printer
        raises(RuntimeError, lambda: runner.run(engine, config), "Unusable folder accepted")
    print("test_printer_errors: OK")


def test_printer_files():
    engine = stochastic.RandomEngine()
    engine.set_seed(3)
    runner = small_builder().build_simulation_runner()
    with tempfile.TemporaryDirectory() as folder:
        printer = stochastic.SimulationPrinterConfigurator(
            output_folder=os.path.join(folder, "sub", "output"),
            statistics_filename="stats.txt",
            statistics_summary_filename="summary.txt",
            realisations_prefix="state_",
            realisations_print_frequency=2,
        )
        config = stochastic.SimulationConfigurator(
            realizations=5, metropolis_hasting_steps=5, burn_in_steps=5, printer=printer
        )
        runner.run(engine, config)
        output = os.path.join(folder, "sub", "output")
        check(
            sorted(os.listdir(output))
            == ["state_0.txt", "state_2.txt", "state_4.txt", "stats.txt", "summary.txt"],
            "Wrong output files: %s" % sorted(os.listdir(output)),
        )
        with open(os.path.join(output, "stats.txt")) as stats:
            lines = stats.read().splitlines()
        # title, names of the terms, one line per realization
        check(len(lines) == 2 + 5, "Wrong statistics lines: %d" % len(lines))
        check("fset_p20" in lines[1], "The terms are not named: " + lines[1])
        # nothing is printed when everything is off
        quiet = stochastic.SimulationPrinterConfigurator(
            output_folder=os.path.join(folder, "quiet"),
            print_statistics=False,
            print_statistics_summary=False,
            print_realisations=False,
        )
        config.printer = quiet
        runner.run(engine, config)
        check(not os.path.exists(os.path.join(folder, "quiet")) or not os.listdir(os.path.join(folder, "quiet")), "Quiet printer wrote files")
        # no printer: the default configurator does not print
        runner.run(engine, stochastic.SimulationConfigurator(realizations=2, metropolis_hasting_steps=1, burn_in_steps=1))
    print("test_printer_files: OK")


def test_printer_field_is_live():
    # printer is a std::optional member: reading it returns a live reference to
    # the stored printer, which cannot be reset to None from Python
    config = stochastic.SimulationConfigurator(
        printer=stochastic.SimulationPrinterConfigurator(
            output_folder="first folder", realisations_print_frequency=7
        )
    )
    printer = config.printer
    check(printer.output_folder == "first folder", "Wrong output folder")
    check(printer.realisations_print_frequency == 7, "Wrong frequency")
    # editing the reference edits the configurator
    printer.output_folder = "changed"
    check(config.printer.output_folder == "changed", "Editing the printer did not change the configurator")
    config.printer.realisations_print_frequency = 3
    check(printer.realisations_print_frequency == 3, "An edit through the field is not seen by the reference")
    # None is refused and the printer is kept
    raises(TypeError, lambda: setattr(config, "printer", None), "printer = None accepted")
    check(config.printer is not None, "The printer was reset")
    check(printer.output_folder == "changed", "The reference changed after a refused printer = None")
    # assigning copies the given printer into the stored one: the reference
    # stays valid and shows the new value
    other = stochastic.SimulationPrinterConfigurator(output_folder="second folder")
    config.printer = other
    check(printer.output_folder == "second folder", "The reference does not show the reassigned printer")
    other.output_folder = "edited copy"
    check(config.printer.output_folder == "second folder", "The assigned printer is shared, not copied")
    # many assignments while holding the reference
    for index in range(100):
        config.printer = stochastic.SimulationPrinterConfigurator(output_folder="folder %d" % index)
        check(printer.output_folder == "folder %d" % index, "The reference lost the printer %d" % index)
    # the reference keeps its configurator alive
    del config
    check(printer.output_folder == "folder 99", "The reference dangles after the configurator is deleted")
    printer.output_folder = "after del"
    check(printer.output_folder == "after del", "The reference cannot be edited after del")
    # a configurator without printer: None by default, None refused, assignable
    empty = stochastic.SimulationConfigurator()
    check(empty.printer is None, "A default configurator has a printer")
    raises(TypeError, lambda: setattr(empty, "printer", None), "printer = None accepted without printer")
    check(empty.printer is None, "A refused None created a printer")
    empty.printer = stochastic.SimulationPrinterConfigurator(output_folder="assigned")
    check(empty.printer.output_folder == "assigned", "Assigning a printer to an empty configurator failed")
    # the keyword constructor still accepts printer=None
    check(stochastic.SimulationConfigurator(printer=None).printer is None, "printer=None refused by the constructor")
    doc = stochastic.SimulationConfigurator.printer.__doc__
    check("live reference" in doc and "None" in doc, "The printer docstring does not say it is live: " + doc)
    check("assign the field" not in doc, "The printer docstring still says it is a copy: " + doc)
    print("test_printer_field_is_live: OK")


if __name__ == "__main__":
    test_printer_files()
    test_printer_errors()
    test_printer_field_is_live()
