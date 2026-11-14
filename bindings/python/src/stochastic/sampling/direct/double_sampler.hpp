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

#pragma once

#include "../../../common.hpp"

#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

namespace geode
{
    void define_double_sampler( pybind11::module &module )
    {
        using Description = DoubleSampler::DistributionDescription;
        pybind11::class_< Description > description( module,
            "DoubleDistributionDescription",
            "Description of a law of a real value: its type and the "
            "parameters it uses (see DoubleSampler.create_distribution). "
            "Parameters that are None are not given." );
        description.def(
            pybind11::init(
                []( std::optional< DistributionType > distribution_type,
                    std::optional< double > min_value,
                    std::optional< double > max_value,
                    std::optional< double > mean,
                    std::optional< double > standard_deviation,
                    std::optional< double > kappa,
                    std::optional< double > alpha,
                    std::optional< std::string > name ) {
                    Description description;
                    detail::set_if_given(
                        description.distribution_type, distribution_type );
                    description.min_value = min_value;
                    description.max_value = max_value;
                    description.mean = mean;
                    description.standard_deviation = standard_deviation;
                    description.kappa = kappa;
                    description.alpha = alpha;
                    detail::set_if_given( description.name, name );
                    return description;
                } ),
            pybind11::arg( "distribution_type" ) = std::nullopt,
            pybind11::kw_only(), pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            pybind11::arg( "mean" ) = std::nullopt,
            pybind11::arg( "standard_deviation" ) = std::nullopt,
            pybind11::arg( "kappa" ) = std::nullopt,
            pybind11::arg( "alpha" ) = std::nullopt,
            pybind11::arg( "name" ) = std::nullopt,
            "Create a description, e.g. DoubleDistributionDescription("
            "\"VonMises\", mean=30.0, kappa=10.0). The default name is "
            "\"default_distribution\" and the default distribution type "
            "is empty (not a valid law)." );
        def_readwrite_named( description, "name", &Description::name,
            "Name of the distribution, used in messages." );
        def_readwrite_named( description, "distribution_type",
            &Description::distribution_type,
            "Type of the law: \"UniformClosed\", \"UniformClosedOpen\", "
            "\"Gaussian\", \"TruncatedGaussian\", \"VonMises\", "
            "\"TruncatedLogNormal\" or \"TruncatedPowerLaw\"." );
        def_readwrite_named( description, "min_value", &Description::min_value,
            "Lower bound (uniform and truncated laws)." );
        def_readwrite_named( description, "max_value", &Description::max_value,
            "Upper bound (uniform and truncated laws)." );
        def_readwrite_named( description, "mean", &Description::mean,
            "Mean (Gaussian, TruncatedGaussian, VonMises, "
            "TruncatedLogNormal)." );
        def_readwrite_named( description, "standard_deviation",
            &Description::standard_deviation,
            "Standard deviation (Gaussian, TruncatedGaussian, "
            "TruncatedLogNormal)." );
        def_readwrite_named( description, "kappa", &Description::kappa,
            "Concentration parameter of the VonMises law." );
        def_readwrite_named( description, "alpha", &Description::alpha,
            "Exponent of the TruncatedPowerLaw law." );
        description
            .def( "string", &Description::string,
                "Return a readable description of the law, one line per "
                "given parameter." )
            .def( "__repr__", []( const Description &self ) {
                return python_call_repr( "DoubleDistributionDescription",
                    { python_repr( self.distribution_type ),
                        "min_value=" + python_repr( self.min_value ),
                        "max_value=" + python_repr( self.max_value ),
                        "mean=" + python_repr( self.mean ),
                        "standard_deviation="
                            + python_repr( self.standard_deviation ),
                        "kappa=" + python_repr( self.kappa ),
                        "alpha=" + python_repr( self.alpha ),
                        "name=" + python_repr( self.name ) } );
            } );

        pybind11::class_< DoubleSampler >( module, "DoubleSampler",
            "Creation of the laws of a real value from their description, and "
            "sampling of these laws." )
            .def_static( "create_distribution",
                &DoubleSampler::create_distribution, pybind11::arg( "desc" ),
                "Create the law described by desc: a UniformClosed, "
                "UniformClosedOpen (min_value and max_value required), "
                "Gaussian (mean and standard_deviation), TruncatedGaussian "
                "(mean and standard_deviation, bounds optional), VonMises "
                "(mean and kappa), TruncatedLogNormal (mean and "
                "standard_deviation, bounds optional) or TruncatedPowerLaw "
                "(alpha, bounds optional). Raise an OpenGeodeDataError if the "
                "type is unknown or a required parameter is missing." )
            .def_static( "create_rad_angle_distribution_from_degree",
                &DoubleSampler::create_rad_angle_distribution_from_degree,
                pybind11::arg( "desc" ),
                "Same as create_distribution for a law of an angle given in "
                "degrees: mean, standard_deviation, min_value and max_value "
                "are converted to radians and normalized to [0, 2*pi) "
                "(kappa and alpha are unchanged)." )
            .def_static( "sample", &DoubleSampler::sample,
                pybind11::arg( "engine" ), pybind11::arg( "dist" ),
                "Sample a value of a law created by create_distribution, "
                "with the engine." );
    }
} // namespace geode
