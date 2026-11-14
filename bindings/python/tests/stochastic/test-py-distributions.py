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

import math

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


def engine_with_seed(seed=7):
    engine = stochastic.RandomEngine()
    engine.set_seed(seed)
    return engine


def test_seed():
    law = stochastic.UniformClosed(min_value=0.0, max_value=10.0)
    first = engine_with_seed(42)
    second = engine_with_seed(42)
    values = [first.sample_uniform_closed(law) for _ in range(20)]
    check(values == [second.sample_uniform_closed(law) for _ in range(20)], "Same seed, different samples")
    other = engine_with_seed(43)
    check(values != [other.sample_uniform_closed(law) for _ in range(20)], "Different seeds, same samples")
    # a string and an integer both seed the engine
    engine_with_seed("a word")
    raises(TypeError, lambda: engine_with_seed(-1), "A negative seed was accepted")
    raises(TypeError, lambda: engine_with_seed(1.5), "A float seed was accepted")
    print("test_seed: OK")


def test_direct_sampling():
    engine = engine_with_seed()
    count = 400
    closed = stochastic.UniformClosed(min_value=2.0, max_value=3.0)
    check(all(2.0 <= engine.sample_uniform_closed(closed) <= 3.0 for _ in range(count)), "UniformClosed range")
    check(engine.sample_uniform_closed(stochastic.UniformClosed(min_value=1.5, max_value=1.5)) == 1.5, "Degenerate range")
    closed_open = stochastic.UniformClosedOpen(min_value=-1.0, max_value=1.0)
    check(all(-1.0 <= engine.sample_uniform_closed_open(closed_open) < 1.0 for _ in range(count)), "UniformClosedOpen range")
    gaussian = stochastic.Gaussian(mean=5.0, standard_deviation=0.1)
    samples = [engine.sample_gaussian(gaussian) for _ in range(count)]
    check(abs(sum(samples) / count - 5.0) < 0.05, "Gaussian mean")
    truncated = stochastic.TruncatedGaussian(mean=0.0, standard_deviation=1.0, min_value=-0.5, max_value=2.0)
    check(all(-0.5 <= engine.sample_truncated_gaussian(truncated) <= 2.0 for _ in range(count)), "TruncatedGaussian range")
    half = stochastic.TruncatedGaussian(min_value=0.0)
    check(all(engine.sample_truncated_gaussian(half) >= 0.0 for _ in range(count)), "Half bounded TruncatedGaussian")
    von_mises = stochastic.VonMises(mean=1.0, concentration=5.0)
    angles = [engine.sample_von_mises(von_mises) for _ in range(count)]
    check(all(0.0 <= angle < 2 * math.pi for angle in angles), "VonMises range")
    check(abs(math.atan2(sum(map(math.sin, angles)), sum(map(math.cos, angles))) - 1.0) < 0.2, "VonMises mean direction")
    lognormal = stochastic.TruncatedLogNormal(mean=0.0, standard_deviation=0.5, min_value=0.5, max_value=3.0)
    check(all(0.5 <= engine.sample_truncated_lognormal(lognormal) <= 3.0 for _ in range(count)), "TruncatedLogNormal range")
    power = stochastic.TruncatedPowerLaw(alpha=2.0, min_value=1.0, max_value=10.0)
    check(all(1.0 <= engine.sample_truncated_powerlaw(power) <= 10.0 for _ in range(count)), "TruncatedPowerLaw range")
    check(all(engine.sample_log() <= 0.0 for _ in range(count)), "sample_log is a logarithm of a probability")
    check(not any(engine.sample_bernoulli(0.0) for _ in range(count)), "Bernoulli(0)")
    check(all(engine.sample_bernoulli(1.0) for _ in range(count)), "Bernoulli(1)")
    mean = sum(engine.sample_bernoulli(0.25) for _ in range(2000)) / 2000
    check(abs(mean - 0.25) < 0.05, "Bernoulli(0.25)")
    print("test_direct_sampling: OK")


def test_sampling_errors():
    engine = engine_with_seed()
    cases = [
        lambda: engine.sample_uniform_closed(stochastic.UniformClosed(min_value=2.0, max_value=1.0)),
        lambda: engine.sample_uniform_closed_open(stochastic.UniformClosedOpen(min_value=1.0, max_value=1.0)),
        lambda: engine.sample_gaussian(stochastic.Gaussian(standard_deviation=0.0)),
        lambda: engine.sample_gaussian(stochastic.Gaussian(mean=float("inf"))),
        lambda: engine.sample_truncated_gaussian(stochastic.TruncatedGaussian(min_value=2.0, max_value=1.0)),
        lambda: engine.sample_von_mises(stochastic.VonMises(concentration=-1.0)),
        lambda: engine.sample_truncated_lognormal(stochastic.TruncatedLogNormal(standard_deviation=-1.0)),
        lambda: engine.sample_truncated_powerlaw(stochastic.TruncatedPowerLaw(alpha=0.0)),
        lambda: engine.sample_bernoulli(1.5),
        lambda: engine.sample_bernoulli(-0.1),
    ]
    for case in cases:
        exception = raises(og.OpenGeodeDataError, case, "No data error for an invalid law")
        check(isinstance(exception, og.OpenGeodeError) and isinstance(exception, RuntimeError), "Error hierarchy")
        check(exception.type == "data", "Error type: " + exception.type)
    print("test_sampling_errors: OK")


def test_is_valid():
    check(stochastic.UniformClosed().is_valid(), "UniformClosed")
    check(stochastic.UniformClosed(min_value=1.0, max_value=1.0).is_valid(), "Equal bounds are valid (warning)")
    check(not stochastic.UniformClosed(min_value=2.0, max_value=1.0).is_valid(), "Inverted bounds")
    check(stochastic.UniformClosedOpen().is_valid(), "UniformClosedOpen")
    check(not stochastic.UniformClosedOpen(min_value=1.0, max_value=1.0).is_valid(), "Empty interval")
    check(stochastic.Gaussian().is_valid(), "Gaussian")
    check(not stochastic.Gaussian(standard_deviation=-1.0).is_valid(), "Negative standard deviation")
    check(not stochastic.Gaussian(mean=float("nan")).is_valid(), "NaN mean")
    check(stochastic.VonMises().is_valid(), "VonMises")
    check(not stochastic.VonMises(concentration=-0.5).is_valid(), "Negative concentration")
    check(stochastic.TruncatedPowerLaw(min_value=1.0, max_value=2.0).is_valid(), "Power law")
    check(not stochastic.TruncatedPowerLaw(alpha=-1.0, min_value=1.0).is_valid(), "Negative exponent")
    check(not stochastic.TruncatedPowerLaw(min_value=3.0, max_value=2.0).is_valid(), "Inverted power law bounds")
    check(stochastic.TruncatedGaussian().is_valid(), "Default truncated gaussian")
    check(stochastic.TruncatedGaussian(mean=10.0, min_value=5.0, max_value=15.0).is_valid(), "Bounded truncated gaussian")
    check(not stochastic.TruncatedGaussian(min_value=7.0).is_valid(), "Truncation 6 sigma above the mean")
    check(not stochastic.TruncatedGaussian(max_value=-7.0).is_valid(), "Truncation 6 sigma below the mean")
    check(not stochastic.TruncatedGaussian(min_value=2.0, max_value=1.0).is_valid(), "Inverted truncation")
    check(not stochastic.TruncatedGaussian(standard_deviation=0.0).is_valid(), "Null truncated gaussian deviation")
    check(not stochastic.TruncatedGaussian(mean=float("inf")).is_valid(), "Infinite truncated gaussian mean")
    check(not stochastic.TruncatedGaussian(min_value=1.0, max_value=1.0).is_valid(), "Empty truncation")
    check(stochastic.TruncatedGaussian(min_value=6.0).is_valid(), "Truncation starting at 6 sigma")
    check(stochastic.TruncatedGaussian(max_value=-6.0).is_valid(), "Truncation ending at -6 sigma")
    check(not stochastic.TruncatedGaussian(max_value=-6.5).is_valid(), "Truncation ending below -6 sigma")
    check(not stochastic.TruncatedLogNormal(min_value=-1.0, max_value=4.0).is_valid(), "Negative log normal lower bound")
    check(stochastic.TruncatedLogNormal(min_value=0.0, max_value=4.0).is_valid(), "Null log normal lower bound")
    check(stochastic.TruncatedLogNormal().is_valid(), "Default log normal")
    check(stochastic.TruncatedLogNormal(min_value=0.5, max_value=3.0).is_valid(), "Bounded log normal")
    check(stochastic.TruncatedLogNormal(min_value=0.5).is_valid(), "Lower bounded log normal")
    check(not stochastic.TruncatedLogNormal(min_value=4.0, max_value=4.0).is_valid(), "Empty log normal range")
    check(not stochastic.TruncatedLogNormal(max_value=0.0).is_valid(), "Log normal upper bound at the default lower bound")
    check(not stochastic.TruncatedLogNormal(standard_deviation=-1.0).is_valid(), "Negative log normal deviation")
    check(not stochastic.TruncatedLogNormal(mean=float("nan")).is_valid(), "NaN log normal mean")
    print("test_is_valid: OK")


def test_double_sampler():
    engine = engine_with_seed()
    description = stochastic.DoubleDistributionDescription(
        "UniformClosed", min_value=1.0, max_value=2.0, name="u"
    )
    law = stochastic.DoubleSampler.create_distribution(description)
    check(isinstance(law, stochastic.UniformClosed), "Wrong law type: " + str(type(law)))
    check((law.min_value, law.max_value) == (1.0, 2.0), "Wrong law parameters")
    check(1.0 <= stochastic.DoubleSampler.sample(engine, law) <= 2.0, "Sample out of range")
    expected = {
        "UniformClosedOpen": (stochastic.UniformClosedOpen, dict(min_value=0.0, max_value=1.0)),
        "Gaussian": (stochastic.Gaussian, dict(mean=0.0, standard_deviation=1.0)),
        "TruncatedGaussian": (stochastic.TruncatedGaussian, dict(mean=0.0, standard_deviation=1.0, max_value=3.0)),
        "VonMises": (stochastic.VonMises, dict(mean=0.5, kappa=2.0)),
        "TruncatedLogNormal": (stochastic.TruncatedLogNormal, dict(mean=0.0, standard_deviation=1.0)),
        "TruncatedPowerLaw": (stochastic.TruncatedPowerLaw, dict(alpha=2.0, min_value=1.0, max_value=5.0)),
    }
    for name, (law_class, parameters) in expected.items():
        law = stochastic.DoubleSampler.create_distribution(
            stochastic.DoubleDistributionDescription(name, **parameters)
        )
        check(isinstance(law, law_class), "%s gave %s" % (name, type(law)))
        value = stochastic.DoubleSampler.sample(engine, law)
        check(isinstance(value, float), name + " sample")
    # parameters of the description reach the law
    law = stochastic.DoubleSampler.create_distribution(
        stochastic.DoubleDistributionDescription("VonMises", mean=0.5, kappa=2.0)
    )
    check((law.mean, law.concentration) == (0.5, 2.0), "VonMises parameters")
    # a law given by the user is accepted too
    check(isinstance(stochastic.DoubleSampler.sample(engine, stochastic.Gaussian()), float), "Sampling a law")
    print("test_double_sampler: OK")


def test_double_sampler_degrees():
    angle = stochastic.DoubleSampler.create_rad_angle_distribution_from_degree(
        stochastic.DoubleDistributionDescription("VonMises", mean=90.0, kappa=10.0)
    )
    check(abs(angle.mean - math.pi / 2) < 1e-9, "Degrees not converted: %r" % angle.mean)
    check(angle.concentration == 10.0, "kappa is not an angle")
    angle = stochastic.DoubleSampler.create_rad_angle_distribution_from_degree(
        stochastic.DoubleDistributionDescription("UniformClosed", min_value=0.0, max_value=180.0)
    )
    check(abs(angle.max_value - math.pi) < 1e-9, "Bounds not converted: %r" % angle.max_value)
    wrapped = stochastic.DoubleSampler.create_rad_angle_distribution_from_degree(
        stochastic.DoubleDistributionDescription("VonMises", mean=450.0, kappa=1.0)
    )
    check(abs(wrapped.mean - math.pi / 2) < 1e-9, "Angle not normalized: %r" % wrapped.mean)
    plain = stochastic.DoubleSampler.create_distribution(
        stochastic.DoubleDistributionDescription("VonMises", mean=90.0, kappa=10.0)
    )
    check(plain.mean == 90.0, "create_distribution must not convert")
    print("test_double_sampler_degrees: OK")


def test_double_sampler_errors():
    create = stochastic.DoubleSampler.create_distribution
    Description = stochastic.DoubleDistributionDescription
    cases = [
        lambda: create(Description("NoSuchLaw", mean=0.0)),
        lambda: create(Description()),
        lambda: create(Description("UniformClosed", min_value=0.0)),
        lambda: create(Description("UniformClosedOpen", max_value=0.0)),
        lambda: create(Description("Gaussian", mean=0.0)),
        lambda: create(Description("TruncatedGaussian", standard_deviation=1.0)),
        lambda: create(Description("VonMises", mean=0.0)),
        lambda: create(Description("TruncatedLogNormal", mean=0.0)),
        lambda: create(Description("TruncatedPowerLaw")),
        lambda: stochastic.DoubleSampler.create_rad_angle_distribution_from_degree(Description("VonMises", kappa=1.0)),
    ]
    for case in cases:
        raises(og.OpenGeodeDataError, case, "No data error for an incomplete description")
    description = Description("Gaussian")
    check(description.name == "default_distribution", "Default name")
    description.name = "renamed"
    check("renamed" in description.string() or description.name == "renamed", "Name setter")
    description.mean = 1.0
    check("mean" in description.string(), "string() lists the given parameters")
    print("test_double_sampler_errors: OK")


if __name__ == "__main__":
    test_seed()
    test_direct_sampling()
    test_sampling_errors()
    test_is_valid()
    test_double_sampler()
    test_double_sampler_degrees()
    test_double_sampler_errors()
