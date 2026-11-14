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

import inspect

import geode_test_env

import opengeode as og
import opengeode_stochastic_py_stochastic as stochastic


def public_classes():
    return [
        (name, value)
        for name, value in vars(stochastic).items()
        if inspect.isclass(value) and not name.startswith("_")
    ]


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def test_documentation_and_argument_names():
    check(stochastic.__doc__, "The module has no docstring")
    classes = public_classes()
    check(len(classes) >= 28, "Classes are missing: %d" % len(classes))
    for class_name, python_class in classes:
        check(python_class.__doc__, class_name + " has no docstring")
        for member_name, member in vars(python_class).items():
            if member_name.startswith("__") and member_name not in (
                "__init__",
                "__eq__",
                "__hash__",
                "__repr__",
            ):
                continue
            documentation = getattr(member, "__doc__", None)
            if isinstance(member, property):
                documentation = member.fget.__doc__
                # the setter signature must name its argument
                if member.fset is not None:
                    check(
                        "arg0" not in (member.fset.__doc__ or "")
                        and "arg1" not in (member.fset.__doc__ or ""),
                        "%s.%s setter has unnamed arguments" % (class_name, member_name),
                    )
            if isinstance(member, staticmethod):
                documentation = member.__func__.__doc__
            if documentation is None:
                continue
            context = class_name + "." + member_name
            check(documentation.strip(), context + " has no docstring")
            check(
                "arg0" not in documentation and "arg1" not in documentation,
                context + " has unnamed arguments:\n" + documentation,
            )
    print("test_documentation_and_argument_names: OK")


def test_every_function_is_documented():
    # a docstring has the signature, then the text
    for class_name, python_class in public_classes():
        for member_name, member in vars(python_class).items():
            if member_name.startswith("_") and member_name != "__init__":
                continue
            if not callable(member) and not isinstance(member, staticmethod):
                continue
            documentation = member.__doc__ or ""
            body = [
                line
                for line in documentation.split("\n")
                if line.strip() and not line.lstrip().startswith(member_name + "(")
                and not line.lstrip().startswith("Overloaded function")
                and not line.lstrip().startswith("1. ")
                and not line.lstrip().startswith("2. ")
            ]
            check(body, "%s.%s has only a signature" % (class_name, member_name))
    print("test_every_function_is_documented: OK")


def test_static_methods():
    for name in ("create_distribution", "create_rad_angle_distribution_from_degree", "sample"):
        check(
            isinstance(inspect.getattr_static(stochastic.DoubleSampler, name), staticmethod),
            "DoubleSampler.%s is not a static method" % name,
        )
    check(
        isinstance(
            inspect.getattr_static(stochastic.OpenGeodeStochasticStochasticLibrary, "initialize"),
            staticmethod,
        ),
        "Library.initialize is not a static method",
    )
    for name in (
        "UniformClosed",
        "UniformClosedOpen",
        "Gaussian",
        "TruncatedGaussian",
        "VonMises",
        "TruncatedLogNormal",
        "TruncatedPowerLaw",
    ):
        law = getattr(stochastic, name)
        check(
            isinstance(inspect.getattr_static(law, "distribution_type_static"), staticmethod),
            name + ".distribution_type_static is not a static method",
        )
        check(
            law.distribution_type_static() == law().distribution_type(),
            name + " static and instance distribution types differ",
        )
        check(law.distribution_type_static().get() == name, name + " has the wrong type name")
    print("test_static_methods: OK")


def test_no_enum_export():
    # no enum in this binding: nothing leaks into the module namespace
    check(not any(name.isupper() for name in dir(stochastic) if not name.startswith("_")), "Exported constants")
    print("test_no_enum_export: OK")


def test_distribution_type_hash_and_equality():
    gaussian = stochastic.DistributionType("Gaussian")
    check(gaussian == stochastic.DistributionType("Gaussian"), "Equal types differ")
    check(hash(gaussian) == hash(stochastic.DistributionType("Gaussian")), "Equal types hash differently")
    check(gaussian != stochastic.DistributionType("VonMises"), "Different types are equal")
    check(gaussian == "Gaussian", "A type is not equal to its name")
    check(gaussian != 3 and gaussian != None, "A type is equal to an unrelated object")  # noqa: E711
    check(gaussian.matches("Gaussian") and not gaussian.matches("VonMises"), "matches is wrong")
    check(len({gaussian, stochastic.DistributionType("Gaussian"), stochastic.DistributionType("")}) == 2, "Set of types")
    check({gaussian: 1}[stochastic.DistributionType("Gaussian")] == 1, "Dict of types")
    check(stochastic.DistributionType().get() == "", "The default type is not empty")
    check(stochastic.Gaussian().distribution_type() == "Gaussian", "Law type")
    print("test_distribution_type_hash_and_equality: OK")


def eval_repr(value):
    namespace = {name: getattr(stochastic, name) for name in dir(stochastic)}
    namespace["og"] = og
    namespace["Point2D"] = og.Point2D
    namespace["Point3D"] = og.Point3D
    namespace["BoundingBox2D"] = og.BoundingBox2D
    namespace["BoundingBox3D"] = og.BoundingBox3D
    return eval(repr(value), namespace)


def test_repr_is_evaluable():
    domain2d = stochastic.SpatialDomainConfig2D(
        min_point=og.Point2D([0.0, 1.0]), max_point=og.Point2D([2.0, 3.0]), buffer_size=0.5
    )
    domain3d = stochastic.SpatialDomainConfig3D(
        min_point=og.Point3D([0.0, 1.0, 2.0]), max_point=og.Point3D([2.0, 3.0, 4.0]), buffer_size=0.25
    )
    description = stochastic.DoubleDistributionDescription(
        "TruncatedGaussian", mean=1.0, standard_deviation=2.0, max_value=float("inf"), name="n"
    )
    sampler = stochastic.Segment2DSamplerConfig(
        length=stochastic.DoubleDistributionDescription("UniformClosed", min_value=0.5, max_value=1.0),
        azimuth=description,
        move_ratio=0.2,
        rotation_ratio=3.0,
    )
    values = [
        stochastic.DistributionType("Gaussian"),
        stochastic.DistributionType(),
        stochastic.UniformClosed(min_value=-1.0, max_value=2.5),
        stochastic.UniformClosedOpen(min_value=-1.0),
        stochastic.Gaussian(mean=1.5, standard_deviation=0.5),
        stochastic.TruncatedGaussian(mean=1.5, min_value=-3.0),
        stochastic.TruncatedGaussian(),
        stochastic.VonMises(mean=0.25, concentration=7.0),
        stochastic.TruncatedLogNormal(standard_deviation=0.3, min_value=1.0, max_value=9.0),
        stochastic.TruncatedPowerLaw(alpha=1.5, min_value=0.5),
        description,
        stochastic.DoubleDistributionDescription(),
        sampler,
        stochastic.ObjectSetDynamicsConfig(birth_ratio=0.5, change_ratio=4.0),
        stochastic.TargetStatisticConfig("a_p20", 12.5, tolerance=0.3),
        stochastic.TargetStatisticConfig(),
        stochastic.SimulationPrinterConfigurator(
            output_folder="/tmp/some folder", print_statistics=False, realisations_print_frequency=7
        ),
        stochastic.SimulationConfigurator(
            realizations=5,
            printer=stochastic.SimulationPrinterConfigurator(output_folder="out"),
        ),
        stochastic.SimulationConfigurator(),
        domain2d,
        domain3d,
    ]
    for value in values:
        text = repr(value)
        check("0x" not in text and "object at" not in text, "Default repr: " + text)
        again = eval_repr(value)
        check(type(again) is type(value), "Wrong type from the repr: " + text)
        check(repr(again) == text, "The repr of %s does not round trip: %s" % (text, repr(again)))
    # nan and inf have no literal
    check("float('inf')" in repr(description), "Infinity repr: " + repr(description))
    check(
        eval_repr(stochastic.Gaussian(mean=float("-inf"))).mean == float("-inf"),
        "Infinite mean lost",
    )
    print("test_repr_is_evaluable: OK")


def test_repr_of_non_value_classes():
    box = og.BoundingBox2D()
    box.add_point(og.Point2D([0.0, 0.0]))
    box.add_point(og.Point2D([4.0, 5.0]))
    domain = stochastic.SpatialDomain2D(box, 1.0)
    # the repr of a bounding box is not evaluable
    check(repr(domain) == "SpatialDomain2D(domain=%r, buffer_size=1.0)" % box, repr(domain))
    builder = stochastic.FractureProcessBuilder()
    check(repr(builder) == "FractureProcessBuilder()", repr(builder))
    check(repr(stochastic.RandomEngine()) == "RandomEngine()", "RandomEngine repr")
    check(repr(stochastic.Segment2DObjectSets()) == "Segment2DObjectSets(nb_sets=0, nb_objects=0)", repr(stochastic.Segment2DObjectSets()))
    check("Point3DObjectSets" in repr(stochastic.Point3DObjectSets()), "Point3DObjectSets repr")
    domain = stochastic.SpatialDomainConfig2D(
        min_point=og.Point2D([0.0, 0.0]), max_point=og.Point2D([5.0, 5.0]), buffer_size=1.0
    )
    builder.set_domain(domain)
    definition = builder.add_fracture_set("fset", 0.1)
    check(repr(definition) == "Segment2DSetDefinition(name='fset', nb_fixed_objects=0)", repr(definition))
    definition.add_fixed_segment(og.Point2D([1.0, 1.0]), og.Point2D([2.0, 2.0]))
    check("nb_fixed_objects=1" in repr(definition), repr(definition))
    print("test_repr_of_non_value_classes: OK")


def test_keyword_constructor_defaults():
    law = stochastic.UniformClosed()
    check((law.min_value, law.max_value) == (0.0, 1.0), "UniformClosed defaults")
    law = stochastic.Gaussian()
    check((law.mean, law.standard_deviation) == (0.0, 1.0), "Gaussian defaults")
    law = stochastic.VonMises()
    check((law.mean, law.concentration) == (0.0, 1.0), "VonMises defaults")
    law = stochastic.TruncatedPowerLaw()
    check((law.alpha, law.min_value, law.max_value) == (2.0, None, None), "Power law defaults")
    law = stochastic.TruncatedGaussian(mean=2.0, max_value=3.0)
    check((law.mean, law.standard_deviation, law.min_value, law.max_value) == (2.0, 1.0, None, 3.0), "Truncated gaussian")
    for positional in (
        lambda: stochastic.Gaussian(1.0, 2.0),
        lambda: stochastic.UniformClosed(0.0, 1.0),
        lambda: stochastic.SimulationConfigurator(5),
    ):
        try:
            positional()
        except TypeError:
            continue
        raise ValueError("[Test] Keyword-only constructor accepted positional arguments")
    # setters validate their type and name their argument
    law = stochastic.Gaussian()
    law.mean = 4.0
    check(law.mean == 4.0, "Setter")
    try:
        law.mean = "not a number"
    except TypeError:
        pass
    else:
        raise ValueError("[Test] Setter accepted a string")
    print("test_keyword_constructor_defaults: OK")


if __name__ == "__main__":
    test_documentation_and_argument_names()
    test_every_function_is_documented()
    test_static_methods()
    test_no_enum_export()
    test_distribution_type_hash_and_equality()
    test_repr_is_evaluable()
    test_repr_of_non_value_classes()
    test_keyword_constructor_defaults()
