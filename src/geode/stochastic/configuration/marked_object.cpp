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

#include <geode/stochastic/configuration/marked_object.hpp>

#include <geode/stochastic/common.hpp>

namespace geode
{
    template < typename Geometry >
    MarkedObject< Geometry >::MarkedObject( Geometry&& geometry )
        : geometry_{ std::move( geometry ) }
    {
    }

    template < typename Geometry >
    MarkedObject< Geometry >::MarkedObject( Geometry&& geometry, Mark mark )
        : geometry_{ std::move( geometry ) }, mark_{ std::move( mark ) }
    {
    }

    template < typename Geometry >
    void MarkedObject< Geometry >::set_geometry( const Geometry& geometry )
    {
        geometry_ = geometry;
    }

    template < typename Geometry >
    void MarkedObject< Geometry >::set_geometry( Geometry&& geometry )
    {
        geometry_ = std::move( geometry );
    }

    template < typename Geometry >
    void MarkedObject< Geometry >::set_mark( const Mark& mark )
    {
        mark_ = mark;
    }

    template < typename Geometry >
    void MarkedObject< Geometry >::set_mark( Mark&& mark )
    {
        mark_ = std::move( mark );
    }

    template < typename Geometry >
    const Geometry& MarkedObject< Geometry >::geometry() const
    {
        return geometry_;
    }

    template < typename Geometry >
    const std::optional< Mark >& MarkedObject< Geometry >::mark() const
    {
        return mark_;
    }

    template < typename Geometry >
    bool MarkedObject< Geometry >::has_mark() const
    {
        return mark_.has_value();
    }

    template class opengeode_stochastic_stochastic_api MarkedObject< Point2D >;
    template class opengeode_stochastic_stochastic_api MarkedObject< Point3D >;
    template class opengeode_stochastic_stochastic_api
        MarkedObject< OwnerSegment2D >;
} // namespace geode