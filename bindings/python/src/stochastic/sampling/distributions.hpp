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

#include <geode/stochastic/sampling/distributions.hpp>

namespace geode
{

    void define_distributions( pybind11::module &module )
    {
        // DistributionType
        pybind11::class_< DistributionType >( module, "DistributionType" )
            .def( pybind11::init<>() )
            .def( pybind11::init< std::string >() )
            .def( "get", &DistributionType::get )
            .def( "matches", &DistributionType::operator== );

        // UniformClosed<double>
        pybind11::class_< UniformClosed< double > >( module, "UniformClosed" )
            .def( pybind11::init<>() )
            .def_readwrite( "min_value", &UniformClosed< double >::min_value )
            .def_readwrite( "max_value", &UniformClosed< double >::max_value )
            .def( "is_valid", &UniformClosed< double >::is_valid )
            .def( "distribution_type_static",
                &UniformClosed< double >::distribution_type_static )
            .def( "distribution_type",
                &UniformClosed< double >::distribution_type )
            .def( "string", &UniformClosed< double >::string );

        // UniformClosedOpen<double>
        pybind11::class_< UniformClosedOpen< double > >(
            module, "UniformClosedOpen" )
            .def( pybind11::init<>() )
            .def_readwrite(
                "min_value", &UniformClosedOpen< double >::min_value )
            .def_readwrite(
                "max_value", &UniformClosedOpen< double >::max_value )
            .def( "is_valid", &UniformClosedOpen< double >::is_valid )
            .def( "distribution_type_static",
                &UniformClosedOpen< double >::distribution_type_static )
            .def( "distribution_type",
                &UniformClosedOpen< double >::distribution_type )
            .def( "string", &UniformClosedOpen< double >::string );

        // Gaussian
        pybind11::class_< Gaussian >( module, "Gaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &Gaussian::mean )
            .def_readwrite(
                "standard_deviation", &Gaussian::standard_deviation )
            .def( "is_valid", &Gaussian::is_valid )
            .def( "distribution_type_static",
                &Gaussian::distribution_type_static )
            .def( "distribution_type", &Gaussian::distribution_type )
            .def( "string", &Gaussian::string );

        // TruncatedGaussian
        pybind11::class_< TruncatedGaussian >( module, "TruncatedGaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &TruncatedGaussian::mean )
            .def_readwrite(
                "standard_deviation", &TruncatedGaussian::standard_deviation )
            .def_readwrite( "min_value", &TruncatedGaussian::min_value )
            .def_readwrite( "max_value", &TruncatedGaussian::max_value )
            .def( "is_valid", &TruncatedGaussian::is_valid )
            .def( "distribution_type_static",
                &TruncatedGaussian::distribution_type_static )
            .def( "distribution_type", &TruncatedGaussian::distribution_type )
            .def( "string", &TruncatedGaussian::string );

        // VonMises
        pybind11::class_< VonMises >( module, "VonMises" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &VonMises::mean )
            .def_readwrite( "concentration", &VonMises::concentration )
            .def( "is_valid", &VonMises::is_valid )
            .def_static( "distribution_type_static",
                &VonMises::distribution_type_static )
            .def( "distribution_type", &VonMises::distribution_type )
            .def( "string", &VonMises::string );

        // TruncatedLogNormal
        pybind11::class_< TruncatedLogNormal >( module, "TruncatedLogNormal" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &TruncatedLogNormal::mean,
                "Mean value of the underlying normal distribution" )
            .def_readwrite( "standard_deviation",
                &TruncatedLogNormal::standard_deviation,
                "Standard deviation value of the underlying normal "
                "distribution" )
            .def_readwrite( "min_value", &TruncatedLogNormal::min_value )
            .def_readwrite( "max_value", &TruncatedLogNormal::max_value )
            .def( "is_valid", &TruncatedLogNormal::is_valid )
            .def_static( "distribution_type_static",
                &TruncatedLogNormal::distribution_type_static )
            .def( "distribution_type", &TruncatedLogNormal::distribution_type )
            .def( "string", &TruncatedLogNormal::string );

        // TruncatedPowerLaw
        pybind11::class_< TruncatedPowerLaw >( module, "TruncatedPowerLaw" )
            .def( pybind11::init<>() )
            .def_readwrite( "alpha", &TruncatedPowerLaw::alpha,
                "Alpha value of the power law" )
            .def_readwrite( "min_value", &TruncatedPowerLaw::min_value )
            .def_readwrite( "max_value", &TruncatedPowerLaw::max_value )
            .def( "is_valid", &TruncatedPowerLaw::is_valid )
            .def_static( "distribution_type_static",
                &TruncatedPowerLaw::distribution_type_static )
            .def( "distribution_type", &TruncatedPowerLaw::distribution_type )
            .def( "string", &TruncatedPowerLaw::string );
    }

} // namespace geode
