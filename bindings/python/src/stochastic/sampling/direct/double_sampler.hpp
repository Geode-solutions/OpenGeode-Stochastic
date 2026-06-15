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

#include <geode/stochastic/sampling/direct/double_sampler.hpp>

namespace geode
{
    void define_double_sampler( pybind11::module &module )
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
            .def_readwrite( "kappa",
                &DoubleSampler::DistributionDescription::kappa,
                "Set up kappa which is the concentration parameter for Von "
                "Mises Distribution law" )
            .def_readwrite( "alpha",
                &DoubleSampler::DistributionDescription::alpha,
                "Set up alpha which is the exponent parameter for power law "
                "Distribution law" )
            .def( "string", &DoubleSampler::DistributionDescription::string,
                "Return a detailed description of the Distribution law" )
            .def( "__repr__",
                []( const DoubleSampler::DistributionDescription &d ) {
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
            .def_static( "create_rad_angle_distribution_from_degree",
                &DoubleSampler::create_rad_angle_distribution_from_degree,
                pybind11::arg( "desc" ),
                "Create a angle distribution in radian from a description "
                "provided in degree" )
            .def_static( "sample", &DoubleSampler::sample,
                pybind11::arg( "engine" ), pybind11::arg( "dist" ),
                "Sample a value from a distribution using a RandomEngine" );
    }
} // namespace geode
