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

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < index_t dimension >
    class SpatialDomain
    {
    public:
        SpatialDomain(
            const BoundingBox< dimension >& domain, double buffer_size )
            : domain_{ domain },
              buffer_size_{ buffer_size },
              extended_domain_{ domain }
        {
            auto volume = domain_.n_volume();
            OpenGeodeStochasticStochasticException::check_exception(
                volume > 0., nullptr, OpenGeodeException::TYPE::data,
                "[SpatialDomain] Undefined Spatial Domain (volume == ", volume,
                ")." );
            OpenGeodeStochasticStochasticException::check_exception(
                buffer_size_ >= 0.0, nullptr, OpenGeodeException::TYPE::data,
                "[SpatialDomain] Buffer size must not be < 0 ( buffer "
                "= ",
                buffer_size_, ")" );
            if( buffer_size_ != 0. )
            {
                extended_domain_.extends( buffer_size_ );
            }
        }

        const BoundingBox< dimension > box() const
        {
            return domain_;
        }

        bool contains( const Point< dimension >& point ) const
        {
            return domain_.contains( point );
        }

        double n_volume() const
        {
            return domain_.n_volume();
        }

        double smallest_length() const
        {
            return std::get< 1 >( domain_.smallest_length() );
        }

        bool extended_contains( const Point< dimension >& point ) const
        {
            return extended_domain_.contains( point );
        }

        double extended_n_volume() const
        {
            return extended_domain_.n_volume();
        }

        const BoundingBox< dimension > extended_box() const
        {
            return extended_domain_;
        }

    private:
        BoundingBox< dimension > domain_;

        double buffer_size_{ 0. };
        BoundingBox< dimension > extended_domain_;
    };

    template < typename ObjectType >
    struct SpatialDomainChecker
    {
        static bool is_anchored_in_domain(
            const SpatialDomain< ObjectType::dim >& domain,
            const ObjectType& obj )
        {
            return domain.contains( obj );
        }

        static bool intersects_domain(
            const SpatialDomain< ObjectType::dim >& domain,
            const ObjectType& obj )
        {
            return domain.contains( obj );
        }
    };

    // Specialization for segment
    template <>
    struct SpatialDomainChecker< OwnerSegment2D >
    {
        static bool is_anchored_in_domain(
            const SpatialDomain< 2 >& domain, const OwnerSegment2D& seg )
        {
            const auto& v = seg.vertices();
            return domain.contains( v[0] );
        }

        static bool intersects_domain(
            const SpatialDomain< 2 >& domain, const OwnerSegment2D& seg )
        {
            return domain.box().intersects( seg );
        }
    };
} // namespace geode