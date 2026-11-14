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

#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < geode::index_t dimension >
    void declare_spatial_domain( pybind11::module &module )
    {
        using Domain = geode::SpatialDomain< dimension >;
        using BBox = geode::BoundingBox< dimension >;

        const auto pyclass_name =
            absl::StrCat( "SpatialDomain", std::to_string( dimension ), "D" );

        pybind11::class_< Domain >( module, pyclass_name.c_str(),
            "Simulation domain: a core domain (the volume of interest), a "
            "buffer zone around it and the extended domain (core domain "
            "extended by the buffer size in every direction)." )
            .def( pybind11::init< const BBox &, double >(),
                pybind11::arg( "domain" ), pybind11::arg( "buffer_size" ),
                "Create a domain from the bounding box of the core domain "
                "and the buffer thickness. Raise an OpenGeodeDataError if "
                "the box has no volume or the buffer size is negative." )
            .def( "box", &Domain::box,
                pybind11::return_value_policy::reference_internal,
                "Return the bounding box of the core domain." )
            .def( "extended_box", &Domain::extended_box,
                pybind11::return_value_policy::reference_internal,
                "Return the bounding box of the extended domain." )
            .def( "buffer_size", &Domain::buffer_size,
                "Return the buffer thickness." )
            .def( "contains", &Domain::contains, pybind11::arg( "point" ),
                "Return True if the point is in the core domain." )
            .def( "extended_contains", &Domain::extended_contains,
                pybind11::arg( "point" ),
                "Return True if the point is in the extended domain." )
            .def( "n_volume", &Domain::n_volume,
                "Return the volume (area in 2D) of the core domain." )
            .def( "extended_n_volume", &Domain::extended_n_volume,
                "Return the volume (area in 2D) of the extended domain." )
            .def( "string", &Domain::string,
                "Return a readable description of the domain." )
            .def( "__repr__", [pyclass_name]( const Domain &self ) {
                return python_call_repr( pyclass_name,
                    { "domain=" + python_repr( self.box() ),
                        "buffer_size=" + python_repr( self.buffer_size() ) } );
            } );
    }

    template < geode::index_t dimension >
    void define_spatial_domain_config( pybind11::module_ &module )
    {
        using DomainConfig = geode::SpatialDomainConfig< dimension >;
        const auto pyclass_name =
            absl::StrCat( "SpatialDomainConfig", dimension, "D" );

        pybind11::class_< DomainConfig > config( module, pyclass_name.c_str(),
            "Description of a box domain with a buffer, to give to the "
            "builder of a simulation." );
        config.def(
            pybind11::init( []( std::optional< Point< dimension > > min_point,
                                std::optional< Point< dimension > > max_point,
                                std::optional< double > buffer_size ) {
                DomainConfig config;
                detail::set_if_given( config.min_point, min_point );
                detail::set_if_given( config.max_point, max_point );
                detail::set_if_given( config.buffer_size, buffer_size );
                return config;
            } ),
            pybind11::kw_only(), pybind11::arg( "min_point" ) = std::nullopt,
            pybind11::arg( "max_point" ) = std::nullopt,
            pybind11::arg( "buffer_size" ) = std::nullopt,
            "Create a configuration (defaults: both corners at the origin, "
            "which is not a valid domain, and buffer_size 0)." );
        def_readwrite_named( config, "min_point", &DomainConfig::min_point,
            "Corner of the box with the smallest coordinates." );
        def_readwrite_named( config, "max_point", &DomainConfig::max_point,
            "Corner of the box with the largest coordinates." );
        def_readwrite_named( config, "buffer_size", &DomainConfig::buffer_size,
            "Thickness of the buffer around the box (>= 0)." );
        config
            .def( "string", &DomainConfig::string,
                "Return a readable description of the configuration." )
            .def( "__repr__", [pyclass_name]( const DomainConfig &self ) {
                return python_call_repr( pyclass_name,
                    { "min_point=" + python_repr( self.min_point ),
                        "max_point=" + python_repr( self.max_point ),
                        "buffer_size=" + python_repr( self.buffer_size ) } );
            } );
    }

    void define_spatial_domain( pybind11::module &module )
    {
        declare_spatial_domain< 2 >( module );
        declare_spatial_domain< 3 >( module );
        define_spatial_domain_config< 2 >( module );
        define_spatial_domain_config< 3 >( module );
    }
} // namespace geode
