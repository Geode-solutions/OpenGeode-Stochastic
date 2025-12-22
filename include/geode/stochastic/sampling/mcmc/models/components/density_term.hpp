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

#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/stochastic/sampling/mcmc/models/components/single_object_term.hpp>

double length_inside_box( const geode::Point2D& p0,
    const geode::Point2D& p1,
    const geode::BoundingBox2D& box )
{
    double dx = p1.value( 0 ) - p0.value( 0 );
    double dy = p1.value( 1 ) - p0.value( 1 );

    double t_min = 0.0;
    double t_max = 1.0;

    auto update_interval = [&]( double p, double q_min, double q_max ) -> bool {
        if( p == 0.0 )
        {
            // Segment parallel to axis
            return ( q_min <= 0.0 && 0.0 <= q_max );
        }

        double t1 = q_min / p;
        double t2 = q_max / p;
        if( t1 > t2 )
            std::swap( t1, t2 );

        t_min = std::max( t_min, t1 );
        t_max = std::min( t_max, t2 );

        return t_min <= t_max;
    };

    // X axis
    if( !update_interval( dx, box.min().value( 0 ) - p0.value( 0 ),
            box.max().value( 0 ) - p0.value( 0 ) ) )
    {
        return 0.0;
    }

    // Y axis
    if( !update_interval( dy, box.min().value( 1 ) - p0.value( 1 ),
            box.max().value( 1 ) - p0.value( 1 ) ) )
    {
        return 0.0;
    }

    if( t_max <= t_min )
    {
        return 0.0;
    }

    double clipped_dx = dx * ( t_max - t_min );
    double clipped_dy = dy * ( t_max - t_min );

    return std::sqrt( clipped_dx * clipped_dx + clipped_dy * clipped_dy );
}

namespace geode
{
    template < typename ObjectType >
    class DensityTerm : public SingleObjectTerm< ObjectType,
                            std::function< double( const ObjectType&,
                                const SpatialDomain< ObjectType::dim >& ) > >
    {
    public:
        explicit DensityTerm( std::string_view name,
            double lambda,
            std::vector< uuid > targeted_set_ids,
            const SpatialDomain< ObjectType::dim >& domain )
            : SingleObjectTerm< ObjectType,
                  std::function< double( const ObjectType&,
                      const SpatialDomain< ObjectType::dim >& ) > >(
                  name,
                  lambda,
                  std::move( targeted_set_ids ),
                  1.0, // scale by domain area to get density per unit
                  []( const ObjectType& obj,
                      const SpatialDomain< ObjectType::dim >& domain ) {
                      if( SpatialDomainChecker<
                              ObjectType >::is_anchored_in_domain( domain,
                              obj ) )
                      {
                          return 1.0;
                      }
                      return 0.0;
                  }, // contribution = 1 anchoredin domain
                  domain )
        {
        }
    };

    class IntensityTerm
        : public SingleObjectTerm< OwnerSegment2D,
              std::function< double( const OwnerSegment2D&,
                  const SpatialDomain< OwnerSegment2D::dim >& ) > >
    {
    public:
        explicit IntensityTerm( std::string_view name,
            double lambda,
            std::vector< uuid > targeted_set_ids,
            double caracteristic_length,
            const SpatialDomain< OwnerSegment2D::dim >& domain )
            : SingleObjectTerm< OwnerSegment2D,
                  std::function< double( const OwnerSegment2D&,
                      const SpatialDomain< OwnerSegment2D::dim >& ) > >(
                  name,
                  lambda,
                  std::move( targeted_set_ids ),
                  1.0 / caracteristic_length,
                  []( const OwnerSegment2D& segment,
                      const SpatialDomain< OwnerSegment2D::dim >& domain ) {
                      auto seg_extremities = segment.vertices();
                      return length_inside_box( seg_extremities[0],
                          seg_extremities[1], domain.box() );
                  },
                  domain )
        {
        }
    };
} // namespace geode