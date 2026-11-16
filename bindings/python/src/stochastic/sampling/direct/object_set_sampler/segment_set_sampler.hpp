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

        /// Property of a distribution description of the sampler config: the
        /// getter gives access to the description in place, the setter
        /// copies the given one and keeps the name of the property.
        template < typename SamplerConfig >
        void def_description_property(
            pybind11::class_< SamplerConfig >& python_class,
            const char* name,
            DoubleSampler::DistributionDescription SamplerConfig::* member,
            const char* doc )
        {
            using Description = DoubleSampler::DistributionDescription;
            const pybind11::cpp_function getter(
                [member]( SamplerConfig& self ) -> Description& {
                    return self.*member;
                },
                pybind11::is_method( python_class ) );
            const pybind11::cpp_function setter(
                [member](
                    SamplerConfig& self, const Description& description ) {
                    assign_keeping_name( self.*member, description );
                },
                pybind11::is_method( python_class ), pybind11::arg( "value" ) );
            python_class.def_property( name, getter, setter,
                pybind11::return_value_policy::reference_internal, doc );
        }
    } // namespace detail

    void define_segment_set_sampler( pybind11::module& module )
    {
        using SamplerConfig = ObjectSamplerConfig< OwnerSegment2D >;
        using Description = DoubleSampler::DistributionDescription;
        pybind11::class_< SamplerConfig > config( module,
            "Segment2DSamplerConfig",
            "Configuration of the stochastic sampler of 2D segments: the "
            "distributions of the length and of the azimuth of a segment "
            "(its first vertex is uniform in the extended domain) and the "
            "relative weights of the moves changing a segment." );
        config.def(
            pybind11::init( []( std::optional< Description > length,
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
            pybind11::arg( "stretch_ratio" ) = std::nullopt,
            "Create a configuration (defaults: move_ratio 0.1, all the "
            "move ratios 1, and length and azimuth distributions that "
            "must be set before the sampler is used). The given "
            "descriptions are copied and keep the names \"length\" and "
            "\"azimuth\"." );
        def_readwrite_named( config, "move_ratio", &SamplerConfig::move_ratio,
            "Relative size of the translation move: the new first vertex is "
            "drawn in a disc of radius move_ratio * segment length." );
        def_readwrite_named( config, "translation_ratio",
            &SamplerConfig::translation_ratio,
            "Relative weight of the translation change move." );
        def_readwrite_named( config, "rotation_ratio",
            &SamplerConfig::rotation_ratio,
            "Relative weight of the change move drawing a new azimuth (first "
            "vertex and length kept)." );
        def_readwrite_named( config, "stretch_ratio",
            &SamplerConfig::stretch_ratio,
            "Relative weight of the change move drawing a new length (first "
            "vertex and azimuth kept)." );
        detail::def_description_property( config, "length",
            &SamplerConfig::length,
            "Distribution of the segment lengths. A max_value (other than for "
            "a Gaussian law) must not exceed the buffer size of the domain, "
            "otherwise building the simulation raises an OpenGeodeDataError; "
            "an unbounded law only logs a warning." );
        detail::def_description_property( config, "azimuth",
            &SamplerConfig::azimuth,
            "Distribution of the segment azimuths, in radians: the direction "
            "of a segment is (sin(azimuth), cos(azimuth)), i.e. the azimuth "
            "is measured clockwise from the y axis. The values are not "
            "converted from degrees." );
        config.def( "__repr__", []( const SamplerConfig& self ) {
            return python_call_repr( "Segment2DSamplerConfig",
                { "length=" + python_repr( self.length ),
                    "azimuth=" + python_repr( self.azimuth ),
                    "move_ratio=" + python_repr( self.move_ratio ),
                    "translation_ratio="
                        + python_repr( self.translation_ratio ),
                    "rotation_ratio=" + python_repr( self.rotation_ratio ),
                    "stretch_ratio=" + python_repr( self.stretch_ratio ) } );
        } );
    }
} // namespace geode
