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

#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace
{
    template < geode::index_t dimension >
    void declare_spatial_domain( pybind11::module &module )
    {
        using Domain = geode::SpatialDomain< dimension >;
        using BBox = geode::BoundingBox< dimension >;

        const auto pyclass_name =
            "SpatialDomain" + std::to_string( dimension ) + "D";

        pybind11::class_< Domain >( module, pyclass_name.c_str() )
            .def( pybind11::init< const BBox &, double >(),
                pybind11::arg( "domain" ), pybind11::arg( "buffer_size" ),
                R"doc(
                Create a spatial domain composed of:

                - a core domain (the VOI)
                - a buffer zone
                - the extended domain (domain + buffer)

                Arguments:
                    domain (BoundingBox): main domain / VOI
                    buffer_size (float): buffer thickness
            )doc" );
    }

    template < geode::index_t dimension >
    void define_spatial_domain_config( pybind11::module_ &module )
    {
        using DomainConfig = geode::SpatialDomainConfig< dimension >;

        pybind11::class_< DomainConfig >( module,
            absl::StrCat( "SpatialDomainConfig", dimension, "D" ).c_str() )
            .def( pybind11::init<>() )
            .def_readwrite( "min_point", &DomainConfig::min_point )
            .def_readwrite( "max_point", &DomainConfig::max_point )
            .def_readwrite( "buffer_size", &DomainConfig::buffer_size )
            .def( "__repr__", []( const DomainConfig &config ) {
                return config.string();
            } );
    }
} // namespace
namespace geode
{
    void define_spatial_domain( pybind11::module &module )
    {
        declare_spatial_domain< 2 >( module );
        declare_spatial_domain< 3 >( module );
        define_spatial_domain_config< 2 >( module );
        define_spatial_domain_config< 3 >( module );
    }
} // namespace geode