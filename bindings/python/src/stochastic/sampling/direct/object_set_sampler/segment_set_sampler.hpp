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

#include "../../../../common.hpp"

#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>

namespace geode
{
    namespace detail
    {
        /// The sampler distributions keep their name ("length", "azimuth")
        /// when they are replaced by a new description.
        inline void assign_keeping_name(
            DoubleSampler::DistributionDescription& field,
            const DoubleSampler::DistributionDescription& description )
        {
            auto name = std::move( field.name );
            field = description;
            field.name = std::move( name );
        }
    } // namespace detail

    void define_segment_set_sampler( pybind11::module& module )
    {
        using SamplerConfig = ObjectSamplerConfig< OwnerSegment2D >;
        using Description = DoubleSampler::DistributionDescription;
        pybind11::class_< SamplerConfig >( module, "Segment2DSamplerConfig",
            "Configuration of the stochastic segment sampler." )
            .def( pybind11::init( []( std::optional< Description > length,
                                      std::optional< Description > azimuth,
                                      std::optional< double > move_ratio,
                                      std::optional< double > translation_ratio,
                                      std::optional< double > rotation_ratio,
                                      std::optional< double > stretch_ratio ) {
                SamplerConfig config;
                if( length )
                {
                    detail::assign_keeping_name( config.length, *length );
                }
                if( azimuth )
                {
                    detail::assign_keeping_name( config.azimuth, *azimuth );
                }
                detail::set_if_given( config.move_ratio, move_ratio );
                detail::set_if_given(
                    config.translation_ratio, translation_ratio );
                detail::set_if_given( config.rotation_ratio, rotation_ratio );
                detail::set_if_given( config.stretch_ratio, stretch_ratio );
                return config;
            } ),
                pybind11::kw_only(), pybind11::arg( "length" ) = std::nullopt,
                pybind11::arg( "azimuth" ) = std::nullopt,
                pybind11::arg( "move_ratio" ) = std::nullopt,
                pybind11::arg( "translation_ratio" ) = std::nullopt,
                pybind11::arg( "rotation_ratio" ) = std::nullopt,
                pybind11::arg( "stretch_ratio" ) = std::nullopt )
            .def_readwrite( "move_ratio", &SamplerConfig::move_ratio,
                "Relative size of the change move (fraction of the segment "
                "length)." )
            .def_readwrite( "translation_ratio",
                &SamplerConfig::translation_ratio,
                "Relative weight of the translation change move." )
            .def_readwrite( "rotation_ratio", &SamplerConfig::rotation_ratio,
                "Relative weight of the change move drawing a new azimuth." )
            .def_readwrite( "stretch_ratio", &SamplerConfig::stretch_ratio,
                "Relative weight of the change move drawing a new length." )
            .def_property(
                "length",
                []( SamplerConfig& self ) -> Description& {
                    return self.length;
                },
                []( SamplerConfig& self, const Description& description ) {
                    detail::assign_keeping_name( self.length, description );
                },
                pybind11::return_value_policy::reference_internal,
                "Distribution used to sample segment lengths." )
            .def_property(
                "azimuth",
                []( SamplerConfig& self ) -> Description& {
                    return self.azimuth;
                },
                []( SamplerConfig& self, const Description& description ) {
                    detail::assign_keeping_name( self.azimuth, description );
                },
                pybind11::return_value_policy::reference_internal,
                "Distribution used to sample segment orientations (degrees, "
                "clockwise from north)." );
    }
} // namespace geode
