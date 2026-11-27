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

#include <geode/stochastic/common.hpp>

#include <geode/basic/range.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < index_t dimension >
    class SpatialDomain
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( SpatialDomain );

    public:
        ~SpatialDomain() = default;

        SpatialDomain( BoundingBox< dimension > domain, double buffer_size )
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

        [[nodiscard]] const BoundingBox< dimension >& box() const
        {
            return domain_;
        }

        [[nodiscard]] bool contains( const Point< dimension >& point ) const
        {
            return domain_.contains( point );
        }

        [[nodiscard]] double n_volume() const
        {
            return domain_.n_volume();
        }

        [[nodiscard]] double smallest_length() const
        {
            return std::get< 1 >( domain_.smallest_length() );
        }

        [[nodiscard]] bool extended_contains(
            const Point< dimension >& point ) const
        {
            return extended_domain_.contains( point );
        }

        [[nodiscard]] double extended_n_volume() const
        {
            return extended_domain_.n_volume();
        }

        [[nodiscard]] double buffer_size() const
        {
            return buffer_size_;
        }

        [[nodiscard]] const BoundingBox< dimension >& extended_box() const
        {
            return extended_domain_;
        }

        [[nodiscard]] std::string string() const
        {
            return absl::StrCat( "Spatial Domain --> center ", domain_.string(),
                " extended: ", extended_domain_.string(),
                " buffer: ", buffer_size_ );
        }

    private:
        BoundingBox< dimension > domain_;

        double buffer_size_{ 0. };
        BoundingBox< dimension > extended_domain_;
    };

    /// Share of an object counted in the domain, the statistic of
    /// ObjectInDomainFeature. A point counts 1 if it is in the closed domain.
    template < typename ObjectType >
    struct SpatialDomainChecker
    {
        [[nodiscard]] static double count_inside_domain(
            const SpatialDomain< ObjectType::dim >& domain,
            const ObjectType& obj )
        {
            return domain.contains( obj ) ? 1. : 0.;
        }

        /// True if the object could be observed in the domain: a point is in
        /// the closed domain.
        [[nodiscard]] static bool touches_domain(
            const SpatialDomain< ObjectType::dim >& domain,
            const ObjectType& obj )
        {
            return domain.contains( obj );
        }
    };

    /// A segment counts 1/2 per endpoint strictly inside the domain, farther
    /// than GLOBAL_EPSILON from its boundary:
    /// - the count does not depend on the endpoint which anchors the segment,
    /// so observed traces need no orientation;
    /// - a trace clipped by the domain ends on its boundary, at a cut point
    /// which is not an endpoint of the fracture: it does not count;
    /// - for a stationary process, the mean count is the mean number of
    /// segments anchored in the domain: lambda * |D|, the density estimator of
    /// Mauldon (1998). A segment crossing the domain without endpoint inside
    /// counts 0: counting the segments hitting the domain would add a term
    /// proportional to its perimeter and to the mean length, which estimates
    /// the mean length, not the density.
    template <>
    struct SpatialDomainChecker< OwnerSegment2D >
    {
        [[nodiscard]] static double count_inside_domain(
            const SpatialDomain< 2 >& domain, const OwnerSegment2D& seg )
        {
            double count{ 0. };
            for( const auto& vertex : seg.vertices() )
            {
                if( is_strictly_inside( domain.box(), vertex ) )
                {
                    count += 0.5;
                }
            }
            return count;
        }

        /// True if the segment could be observed in the domain: it intersects
        /// the closed domain, wherever its endpoints are.
        [[nodiscard]] static bool touches_domain(
            const SpatialDomain< 2 >& domain, const OwnerSegment2D& seg )
        {
            return domain.box().intersects( seg );
        }

    private:
        [[nodiscard]] static bool is_strictly_inside(
            const BoundingBox2D& box, const Point2D& point )
        {
            for( const auto axis : LRange{ 2 } )
            {
                if( point.value( axis )
                        <= box.min().value( axis ) + GLOBAL_EPSILON
                    || point.value( axis )
                           >= box.max().value( axis ) - GLOBAL_EPSILON )
                {
                    return false;
                }
            }
            return true;
        }
    };

    template < index_t dimension >
    struct SpatialDomainConfig
    {
        Point< dimension > min_point;
        Point< dimension > max_point;
        double buffer_size{ 0.0 };

        [[nodiscard]] std::string string() const
        {
            std::string message = "Domain: ";
            absl::StrAppend(
                &message, "\n\t --> Min point: ", min_point.string() );
            absl::StrAppend(
                &message, "\n\t --> Max point: ", max_point.string() );
            absl::StrAppend( &message, "\n\t --> Buffer size: ", buffer_size );
            return message;
        }
    };

    template < index_t dimension >
    std::unique_ptr< SpatialDomain< dimension > > build_spatial_domain(
        const SpatialDomainConfig< dimension >& config )
    {
        BoundingBox< dimension > box{ config.min_point, config.max_point };

        return std::make_unique< SpatialDomain< dimension > >(
            std::move( box ), config.buffer_size );
    }
} // namespace geode