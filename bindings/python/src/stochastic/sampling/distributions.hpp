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

#include "../../common.hpp"

#include <geode/stochastic/sampling/distributions.hpp>

namespace geode
{
    namespace detail
    {
        /// Static description of a distribution, repeated by the classes
        template < typename Law >
        void define_distribution_common(
            pybind11::class_< Law > &python_class, const char *doc_is_valid )
        {
            python_class.def( "is_valid", &Law::is_valid, doc_is_valid )
                .def_static( "distribution_type_static",
                    &Law::distribution_type_static,
                    "Return the DistributionType shared by all the laws of "
                    "this class." )
                .def( "distribution_type", &Law::distribution_type,
                    "Return the DistributionType of the law." )
                .def( "string", &Law::string,
                    "Return a readable description of the law." );
        }
    } // namespace detail

    void define_distributions( pybind11::module &module )
    {
        pybind11::class_< DistributionType >( module, "DistributionType",
            "Name of a distribution law, e.g. \"Gaussian\" or "
            "\"VonMises\". A str is accepted wherever a DistributionType "
            "is expected." )
            .def( pybind11::init<>(),
                "Create an empty distribution type (no law)." )
            .def( pybind11::init< std::string >(), pybind11::arg( "name" ),
                "Create the distribution type of the given name." )
            .def( "get", &DistributionType::get,
                "Return the name of the distribution type." )
            .def( "matches", &DistributionType::operator==,
                pybind11::arg( "other" ),
                "Return True if other has the same name." )
            .def( "__eq__", &DistributionType::operator==,
                pybind11::is_operator(), pybind11::arg( "other" ),
                "Return True if other has the same name." )
            .def(
                "__hash__",
                []( const DistributionType &self ) {
                    return std::hash< std::string >{}( self.get() );
                },
                "Return a hash of the name." )
            .def( "__repr__", []( const DistributionType &self ) {
                return python_call_repr(
                    "DistributionType", { python_repr( self.get() ) } );
            } );
        // a distribution type can be given as a string
        pybind11::implicitly_convertible< std::string, DistributionType >();

        using UniformClosedLaw = UniformClosed< double >;
        pybind11::class_< UniformClosedLaw > uniform_closed( module,
            "UniformClosed",
            "Uniform law on the closed interval [min_value, max_value]; "
            "sampling requires min_value <= max_value." );
        uniform_closed.def(
            pybind11::init( []( std::optional< double > min_value,
                                std::optional< double > max_value ) {
                UniformClosedLaw law;
                detail::set_if_given( law.min_value, min_value );
                detail::set_if_given( law.max_value, max_value );
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            "Create the law (defaults: min_value 0, max_value 1)." );
        def_readwrite_named( uniform_closed, "min_value",
            &UniformClosedLaw::min_value, "Lower bound of the interval." );
        def_readwrite_named( uniform_closed, "max_value",
            &UniformClosedLaw::max_value, "Upper bound of the interval." );
        detail::define_distribution_common( uniform_closed,
            "Return True if min_value <= max_value (a warning is logged when "
            "they are equal, an error when min_value > max_value)." );
        uniform_closed.def( "__repr__", []( const UniformClosedLaw &law ) {
            return python_call_repr( "UniformClosed",
                { "min_value=" + python_repr( law.min_value ),
                    "max_value=" + python_repr( law.max_value ) } );
        } );

        using UniformClosedOpenLaw = UniformClosedOpen< double >;
        pybind11::class_< UniformClosedOpenLaw > uniform_closed_open( module,
            "UniformClosedOpen",
            "Uniform law on the half-open interval [min_value, max_value); "
            "sampling requires min_value < max_value." );
        uniform_closed_open.def(
            pybind11::init( []( std::optional< double > min_value,
                                std::optional< double > max_value ) {
                UniformClosedOpenLaw law;
                detail::set_if_given( law.min_value, min_value );
                detail::set_if_given( law.max_value, max_value );
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            "Create the law (defaults: min_value 0, max_value 1)." );
        def_readwrite_named( uniform_closed_open, "min_value",
            &UniformClosedOpenLaw::min_value,
            "Lower bound of the interval (included)." );
        def_readwrite_named( uniform_closed_open, "max_value",
            &UniformClosedOpenLaw::max_value,
            "Upper bound of the interval (excluded)." );
        detail::define_distribution_common( uniform_closed_open,
            "Return True if min_value < max_value (an error is logged "
            "otherwise)." );
        uniform_closed_open.def(
            "__repr__", []( const UniformClosedOpenLaw &law ) {
                return python_call_repr( "UniformClosedOpen",
                    { "min_value=" + python_repr( law.min_value ),
                        "max_value=" + python_repr( law.max_value ) } );
            } );

        pybind11::class_< Gaussian > gaussian( module, "Gaussian",
            "Gaussian law; sampling requires a finite mean and a finite "
            "standard_deviation > 0." );
        gaussian.def(
            pybind11::init( []( std::optional< double > mean,
                                std::optional< double > standard_deviation ) {
                Gaussian law;
                detail::set_if_given( law.mean, mean );
                detail::set_if_given(
                    law.standard_deviation, standard_deviation );
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "mean" ) = std::nullopt,
            pybind11::arg( "standard_deviation" ) = std::nullopt,
            "Create the law (defaults: mean 0, standard_deviation 1)." );
        def_readwrite_named( gaussian, "mean", &Gaussian::mean, "Mean." );
        def_readwrite_named( gaussian, "standard_deviation",
            &Gaussian::standard_deviation, "Standard deviation." );
        detail::define_distribution_common( gaussian,
            "Return True if the mean is finite and the standard deviation is "
            "finite and > 0 (an error is logged otherwise)." );
        gaussian.def( "__repr__", []( const Gaussian &law ) {
            return python_call_repr(
                "Gaussian", { "mean=" + python_repr( law.mean ),
                                "standard_deviation="
                                    + python_repr( law.standard_deviation ) } );
        } );

        pybind11::class_< TruncatedGaussian > truncated_gaussian( module,
            "TruncatedGaussian",
            "Gaussian law restricted to [min_value, max_value]; a bound that "
            "is None is infinite. Sampling requires a finite mean, a finite "
            "standard_deviation > 0 and min_value < max_value." );
        truncated_gaussian.def(
            pybind11::init( []( std::optional< double > mean,
                                std::optional< double > standard_deviation,
                                std::optional< double > min_value,
                                std::optional< double > max_value ) {
                TruncatedGaussian law;
                detail::set_if_given( law.mean, mean );
                detail::set_if_given(
                    law.standard_deviation, standard_deviation );
                law.min_value = min_value;
                law.max_value = max_value;
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "mean" ) = std::nullopt,
            pybind11::arg( "standard_deviation" ) = std::nullopt,
            pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            "Create the law (defaults: mean 0, standard_deviation 1, no "
            "bounds)." );
        def_readwrite_named( truncated_gaussian, "mean",
            &TruncatedGaussian::mean, "Mean of the untruncated law." );
        def_readwrite_named( truncated_gaussian, "standard_deviation",
            &TruncatedGaussian::standard_deviation,
            "Standard deviation of the untruncated law." );
        def_readwrite_named( truncated_gaussian, "min_value",
            &TruncatedGaussian::min_value, "Lower bound, None for -infinity." );
        def_readwrite_named( truncated_gaussian, "max_value",
            &TruncatedGaussian::max_value, "Upper bound, None for +infinity." );
        detail::define_distribution_common( truncated_gaussian,
            "Return True if the parameters pass the checks of the C++ "
            "TruncatedGaussian::is_valid (an error is logged otherwise)." );
        truncated_gaussian.def( "__repr__", []( const TruncatedGaussian &law ) {
            return python_call_repr( "TruncatedGaussian",
                { "mean=" + python_repr( law.mean ),
                    "standard_deviation="
                        + python_repr( law.standard_deviation ),
                    "min_value=" + python_repr( law.min_value ),
                    "max_value=" + python_repr( law.max_value ) } );
        } );

        pybind11::class_< VonMises > von_mises( module, "VonMises",
            "Von Mises law of an angle in radians: mean direction and "
            "concentration kappa. Sampling requires a finite mean and a "
            "concentration >= 0, and gives an angle in [0, 2*pi)." );
        von_mises.def(
            pybind11::init( []( std::optional< double > mean,
                                std::optional< double > concentration ) {
                VonMises law;
                detail::set_if_given( law.mean, mean );
                detail::set_if_given( law.concentration, concentration );
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "mean" ) = std::nullopt,
            pybind11::arg( "concentration" ) = std::nullopt,
            "Create the law (defaults: mean 0, concentration 1)." );
        def_readwrite_named(
            von_mises, "mean", &VonMises::mean, "Mean direction in radians." );
        def_readwrite_named( von_mises, "concentration",
            &VonMises::concentration, "Concentration parameter kappa (>= 0)." );
        detail::define_distribution_common( von_mises,
            "Return True if the concentration is >= 0 and the mean is "
            "finite." );
        von_mises.def( "__repr__", []( const VonMises &law ) {
            return python_call_repr( "VonMises",
                { "mean=" + python_repr( law.mean ),
                    "concentration=" + python_repr( law.concentration ) } );
        } );

        pybind11::class_< TruncatedLogNormal > truncated_lognormal( module,
            "TruncatedLogNormal",
            "Log-normal law restricted to [min_value, max_value]; mean and "
            "standard_deviation are the ones of the underlying normal law. "
            "Sampling requires a finite mean, a finite standard_deviation > 0 "
            "and min_value < max_value (None means 0 for the lower bound and "
            "+infinity for the upper bound)." );
        truncated_lognormal.def(
            pybind11::init( []( std::optional< double > mean,
                                std::optional< double > standard_deviation,
                                std::optional< double > min_value,
                                std::optional< double > max_value ) {
                TruncatedLogNormal law;
                detail::set_if_given( law.mean, mean );
                detail::set_if_given(
                    law.standard_deviation, standard_deviation );
                law.min_value = min_value;
                law.max_value = max_value;
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "mean" ) = std::nullopt,
            pybind11::arg( "standard_deviation" ) = std::nullopt,
            pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            "Create the law (defaults: mean 0, standard_deviation 1, no "
            "bounds)." );
        def_readwrite_named( truncated_lognormal, "mean",
            &TruncatedLogNormal::mean, "Mean of the underlying normal law." );
        def_readwrite_named( truncated_lognormal, "standard_deviation",
            &TruncatedLogNormal::standard_deviation,
            "Standard deviation of the underlying normal law." );
        def_readwrite_named( truncated_lognormal, "min_value",
            &TruncatedLogNormal::min_value, "Lower bound, None for 0." );
        def_readwrite_named( truncated_lognormal, "max_value",
            &TruncatedLogNormal::max_value,
            "Upper bound, None for +infinity." );
        detail::define_distribution_common( truncated_lognormal,
            "Return True if the parameters pass the checks of the C++ "
            "TruncatedLogNormal::is_valid (an error is logged otherwise)." );
        truncated_lognormal.def(
            "__repr__", []( const TruncatedLogNormal &law ) {
                return python_call_repr( "TruncatedLogNormal",
                    { "mean=" + python_repr( law.mean ),
                        "standard_deviation="
                            + python_repr( law.standard_deviation ),
                        "min_value=" + python_repr( law.min_value ),
                        "max_value=" + python_repr( law.max_value ) } );
            } );

        pybind11::class_< TruncatedPowerLaw > truncated_powerlaw( module,
            "TruncatedPowerLaw",
            "Power law of exponent alpha restricted to [min_value, "
            "max_value]. Sampling requires alpha > 0 and min_value < "
            "max_value (None means a tiny positive lower bound and "
            "+infinity for the upper bound)." );
        truncated_powerlaw.def(
            pybind11::init( []( std::optional< double > alpha,
                                std::optional< double > min_value,
                                std::optional< double > max_value ) {
                TruncatedPowerLaw law;
                detail::set_if_given( law.alpha, alpha );
                law.min_value = min_value;
                law.max_value = max_value;
                return law;
            } ),
            pybind11::kw_only(), pybind11::arg( "alpha" ) = std::nullopt,
            pybind11::arg( "min_value" ) = std::nullopt,
            pybind11::arg( "max_value" ) = std::nullopt,
            "Create the law (defaults: alpha 2, no bounds)." );
        def_readwrite_named( truncated_powerlaw, "alpha",
            &TruncatedPowerLaw::alpha, "Exponent of the power law (> 0)." );
        def_readwrite_named( truncated_powerlaw, "min_value",
            &TruncatedPowerLaw::min_value,
            "Lower bound, None for a tiny positive value." );
        def_readwrite_named( truncated_powerlaw, "max_value",
            &TruncatedPowerLaw::max_value, "Upper bound, None for +infinity." );
        detail::define_distribution_common( truncated_powerlaw,
            "Return True if alpha is finite and > 0 and the bounds are "
            "positive and ordered (an error is logged otherwise)." );
        truncated_powerlaw.def( "__repr__", []( const TruncatedPowerLaw &law ) {
            return python_call_repr( "TruncatedPowerLaw",
                { "alpha=" + python_repr( law.alpha ),
                    "min_value=" + python_repr( law.min_value ),
                    "max_value=" + python_repr( law.max_value ) } );
        } );
    }

} // namespace geode
