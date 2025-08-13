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

#pragma once

#include <optional>

#include <geode/stochastic/common.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( Point );
} // namespace geode
namespace geode
{
    struct Mark
    {
        int label;

        bool operator==( const Mark& other ) const
        {
            return label == other.label;
        }
    };

    template < typename Geometry >
    class MarkedObject
    {
    public:
        MarkedObject( Geometry&& geometry );
        MarkedObject( Geometry&& geometry, Mark mark );

        void set_geometry( const Geometry& geometry );
        void set_geometry( Geometry&& geometry );

        void set_mark( const Mark& mark );
        void set_mark( Mark&& mark );

        const Geometry& geometry() const;
        const std::optional< Mark >& mark() const;
        bool has_mark() const;

        auto bounding_box() const
        {
            if constexpr( std::is_same_v< Geometry, Point2D > )
            {
                geode::BoundingBox< 2 > box;
                box.add_point( geometry_ );
                return box;
            }
            else if constexpr( std::is_same_v< Geometry, Point3D > )
            {
                geode::BoundingBox< 3 > box;
                box.add_point( geometry_ );
                return box;
            }
            else
            {
                return geometry_.bounding_box();
            }
        }

        auto barycenter() const
        {
            if constexpr( std::is_same_v< Geometry,
                              Point2D > || std::is_same_v< Geometry, Point3D > )
            {
                return geometry_;
            }
            else
            {
                return geometry_.barycenter();
            }
        }

    private:
        Geometry geometry_;
        std::optional< Mark > mark_;
    };
} // namespace geode