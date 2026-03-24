/*
 * Copyright (c) 2019 - 2026 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <geode/stochastic/sampling/mcmc/energy_terms/intensity_term.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>

namespace
{
    double length_inside_box( const geode::Point2D& segment_start,
        const geode::Point2D& segment_end,
        const geode::BoundingBox2D& box )
    {
        // Segment direction
        const double dir_x = segment_end.value( 0 ) - segment_start.value( 0 );
        const double dir_y = segment_end.value( 1 ) - segment_start.value( 1 );

        // Parameter interval where the segment is inside the box
        double t_enter = 0.0; // start of valid interval
        double t_exit = 1.0; // end of valid interval

        // Clips the parametric segment against one axis interval
        auto clip_against_axis_interval = [&t_enter, &t_exit]( double direction,
                                              double min_distance,
                                              double max_distance ) -> bool {
            // Segment is parallel to this axis
            if( std::fabs( direction ) < geode::GLOBAL_EPSILON )
            {
                // Outside the interval → no intersection
                return min_distance <= 0.0 && 0.0 <= max_distance;
            }

            double axis_min_t = min_distance / direction;
            double axis_max_t = max_distance / direction;

            if( axis_min_t > axis_max_t )
            {
                std::swap( axis_min_t, axis_max_t );
            }

            t_enter = std::max( t_enter, axis_min_t );
            t_exit = std::min( t_exit, axis_max_t );

            return t_enter <= t_exit;
        };

        // Clip against X interval of the box
        if( !clip_against_axis_interval( dir_x,
                box.min().value( 0 ) - segment_start.value( 0 ),
                box.max().value( 0 ) - segment_start.value( 0 ) ) )
        {
            return 0.0;
        }

        // Clip against Y interval of the box
        if( !clip_against_axis_interval( dir_y,
                box.min().value( 1 ) - segment_start.value( 1 ),
                box.max().value( 1 ) - segment_start.value( 1 ) ) )
        {
            return 0.0;
        }

        // No portion of the segment is inside the box
        if( t_exit <= t_enter )
        {
            return 0.0;
        }

        // Length of the clipped segment
        const double clipped_dx = dir_x * ( t_exit - t_enter );
        const double clipped_dy = dir_y * ( t_exit - t_enter );

        return std::sqrt(
            ( clipped_dx * clipped_dx ) + ( clipped_dy * clipped_dy ) );
    }
} // namespace
namespace geode
{
    IntensityTerm::IntensityTerm( std::string_view name,
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
                  const SpatialDomain< OwnerSegment2D::dim >& spatial_domain ) {
                  auto seg_extremities = segment.vertices();
                  return length_inside_box( seg_extremities[0],
                      seg_extremities[1], spatial_domain.box() );
              },
              domain )
    {
    }
} // namespace geode