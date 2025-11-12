/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/stochastic/sampling/direct/double_sampler.hpp>

namespace geode
{
    void define_double_sampler( pybind11::module& module )
    {
        // Bind DistributionDescription
        pybind11::class_< DoubleSampler::DistributionDescription >(
            module, "DoubleDistributionDescription" )
            .def( pybind11::init<>() )
            .def_readwrite(
                "name", &DoubleSampler::DistributionDescription::name )
            .def_readwrite( "distribution_type",
                &DoubleSampler::DistributionDescription::distribution_type )
            .def_readwrite( "min_value",
                &DoubleSampler::DistributionDescription::min_value )
            .def_readwrite( "max_value",
                &DoubleSampler::DistributionDescription::max_value )
            .def_readwrite(
                "mean", &DoubleSampler::DistributionDescription::mean )
            .def_readwrite( "standard_deviation",
                &DoubleSampler::DistributionDescription::standard_deviation )
            .def( "string", &DoubleSampler::DistributionDescription::string,
                "Return a detailed description of the Distribution law" )
            .def( "__repr__",
                []( const DoubleSampler::DistributionDescription& d ) {
                    return "<DoubleDistributionDescription name=" + d.name
                           + ">";
                } );

        // Bind Distribution variant
        pybind11::class_< DoubleSampler::Distribution >(
            module, "Distribution" );

        // Bind DoubleSampler
        pybind11::class_< DoubleSampler >( module, "DoubleSampler" )
            .def_static( "create_distribution",
                &DoubleSampler::create_distribution, pybind11::arg( "desc" ),
                "Create a distribution from a description" )
            .def_static( "sample", &DoubleSampler::sample,
                pybind11::arg( "engine" ), pybind11::arg( "dist" ),
                "Sample a value from a distribution using a RandomEngine" );

        // Optionally, expose the variant types
        pybind11::class_< UniformClosed< double > >(
            module, "UniformClosedDouble" )
            .def( pybind11::init<>() )
            .def_readwrite( "min_value", &UniformClosed< double >::min_value )
            .def_readwrite( "max_value", &UniformClosed< double >::max_value );

        pybind11::class_< UniformClosedOpen< double > >(
            module, "UniformClosedOpenDouble" )
            .def( pybind11::init<>() )
            .def_readwrite(
                "min_value", &UniformClosedOpen< double >::min_value )
            .def_readwrite(
                "max_value", &UniformClosedOpen< double >::max_value );

        pybind11::class_< Gaussian >( module, "Gaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &Gaussian::mean )
            .def_readwrite(
                "standard_deviation", &Gaussian::standard_deviation );

        pybind11::class_< TruncatedGaussian >( module, "TruncatedGaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &TruncatedGaussian::mean )
            .def_readwrite(
                "standard_deviation", &TruncatedGaussian::standard_deviation )
            .def_readwrite( "min_value", &TruncatedGaussian::min_value )
            .def_readwrite( "max_value", &TruncatedGaussian::max_value );
    }
} // namespace geode
