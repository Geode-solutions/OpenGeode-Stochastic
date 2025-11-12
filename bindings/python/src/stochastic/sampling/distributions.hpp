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

#include <geode/stochastic/sampling/distributions.hpp>

namespace geode
{

    void define_distributions( pybind11::module& module )
    {
        // DistributionType
        pybind11::class_< DistributionType >( module, "DistributionType" )
            .def( pybind11::init<>() )
            .def( pybind11::init< std::string >() )
            .def( "get", &DistributionType::get )
            .def( "matches", &DistributionType::operator== );

        // UniformClosed<double>
        pybind11::class_< geode::UniformClosed< double > >( m, "UniformClosed" )
            .def( pybind11::init<>() )
            .def_readwrite(
                "min_value", &geode::UniformClosed< double >::min_value )
            .def_readwrite(
                "max_value", &geode::UniformClosed< double >::max_value )
            .def( "is_valid", &geode::UniformClosed< double >::is_valid )
            .def( "distribution_type_static",
                &geode::UniformClosed< double >::distribution_type_static )
            .def( "distribution_type",
                &geode::UniformClosed< double >::distribution_type )
            .def( "string", &geode::UniformClosed< double >::string );

        // UniformClosedOpen<double>
        pybind11::class_< geode::UniformClosedOpen< double > >(
            m, "UniformClosedOpen" )
            .def( pybind11::init<>() )
            .def_readwrite(
                "min_value", &geode::UniformClosedOpen< double >::min_value )
            .def_readwrite(
                "max_value", &geode::UniformClosedOpen< double >::max_value )
            .def( "is_valid", &geode::UniformClosedOpen< double >::is_valid )
            .def( "distribution_type_static",
                &geode::UniformClosedOpen< double >::distribution_type_static )
            .def( "distribution_type",
                &geode::UniformClosedOpen< double >::distribution_type )
            .def( "string", &geode::UniformClosedOpen< double >::string );

        // Gaussian
        pybind11::class_< geode::Gaussian >( m, "Gaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &geode::Gaussian::mean )
            .def_readwrite(
                "standard_deviation", &geode::Gaussian::standard_deviation )
            .def( "is_valid", &geode::Gaussian::is_valid )
            .def( "distribution_type_static",
                &geode::Gaussian::distribution_type_static )
            .def( "distribution_type", &geode::Gaussian::distribution_type )
            .def( "string", &geode::Gaussian::string );

        // TruncatedGaussian
        pybind11::class_< geode::TruncatedGaussian >( m, "TruncatedGaussian" )
            .def( pybind11::init<>() )
            .def_readwrite( "mean", &geode::TruncatedGaussian::mean )
            .def_readwrite( "standard_deviation",
                &geode::TruncatedGaussian::standard_deviation )
            .def_readwrite( "min_value", &geode::TruncatedGaussian::min_value )
            .def_readwrite( "max_value", &geode::TruncatedGaussian::max_value )
            .def( "is_valid", &geode::TruncatedGaussian::is_valid )
            .def( "distribution_type_static",
                &geode::TruncatedGaussian::distribution_type_static )
            .def( "distribution_type",
                &geode::TruncatedGaussian::distribution_type )
            .def( "string", &geode::TruncatedGaussian::string );
    }

} // namespace geode
