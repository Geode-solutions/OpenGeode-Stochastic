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
#include <geode/stochastic/common.hpp>

#include <geode/basic/assert.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

geode::SpatialDomain< 2 > init_domain()
{
    geode::BoundingBox2D box;
    box.add_point( geode::Point2D{ { 0., 0. } } );
    box.add_point( geode::Point2D{ { 1.0, 1.0 } } );
    return geode::SpatialDomain< 2 >{ box, 0.5 };
}

void test_spatial_domain_2D()
{
    using namespace geode;

    auto domain = init_domain();

    // Test volumes
    OPENGEODE_EXCEPTION(
        domain.n_volume() == 1.0, "[Test] Domain volume wrong." );
    OPENGEODE_EXCEPTION( domain.extended_n_volume() == 4,
        "[Test] Extended domain volume wrong." ); // (1+0.5*2)^2

    // Test points inside domain
    geode::Point2D inside{ { 0.5, 0.5 } };
    geode::Point2D boundary{ { 1., 1. } };
    geode::Point2D outside{ { 1.4, 1.4 } };

    OPENGEODE_EXCEPTION(
        domain.contains( inside ), "[Test] Point inside should be contained." );
    OPENGEODE_EXCEPTION( domain.contains( boundary ),
        "[Test] Boundary point should be contained." );
    OPENGEODE_EXCEPTION( !domain.contains( outside ),
        "[Test] Point outside should not be contained." );

    OPENGEODE_EXCEPTION( domain.extended_contains( outside ),
        "[Test] Point outside original but in buffer should be contained." );

    // Test SpatialDomainChecker for points
    OPENGEODE_EXCEPTION( SpatialDomainChecker< Point2D >::is_anchored_in_domain(
                             domain, inside ),
        "[Test] anchored_in_domain failed for inside point" );
    OPENGEODE_EXCEPTION(
        !SpatialDomainChecker< Point2D >::is_anchored_in_domain(
            domain, outside ),
        "[Test] anchored_in_domain failed for outside point" );
    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< Point2D >::intersects_domain( domain, boundary ),
        "[Test] intersects_domain failed for boundary point" );

    // Test segments
    geode::OwnerSegment2D seg_inside{ geode::Point2D{ { 0.1, 0.1 } },
        geode::Point2D{ { 0.9, 0.9 } } };
    geode::OwnerSegment2D seg_partial_anchored{ geode::Point2D{ { 0.9, 0.9 } },
        geode::Point2D{ { 1.6, 1.6 } } };
    geode::OwnerSegment2D seg_partial{ geode::Point2D{ { -0.5, 0.9 } },
        geode::Point2D{ { 0.6, 0.6 } } };
    geode::OwnerSegment2D seg_outside{ geode::Point2D{ { 1.7, 1.7 } },
        geode::Point2D{ { 2.0, 2.0 } } };
    geode::OwnerSegment2D seg_cross{ geode::Point2D{ { -0.5, 0.5 } },
        geode::Point2D{ { 1.5, 0.5 } } };

    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< OwnerSegment2D >::is_anchored_in_domain(
            domain, seg_inside ),
        "[Test] Segment inside VOI should be anchored." );
    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< OwnerSegment2D >::is_anchored_in_domain(
            domain, seg_partial_anchored ),
        "[Test] Segment with lower left extremity inside VOI should be "
        "anchored." );
    OPENGEODE_EXCEPTION(
        !SpatialDomainChecker< OwnerSegment2D >::is_anchored_in_domain(
            domain, seg_partial ),
        "[Test] Segment with lower left extremity outside VOI is not anchored "
        "by definition." );
    OPENGEODE_EXCEPTION(
        !SpatialDomainChecker< OwnerSegment2D >::is_anchored_in_domain(
            domain, seg_outside ),
        "[Test] Segment completely outside should not be anchored." );

    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< OwnerSegment2D >::intersects_domain(
            domain, seg_inside ),
        "[Test] Segment inside VOI should intersect." );
    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< OwnerSegment2D >::intersects_domain(
            domain, seg_partial ),
        "[Test] Segment partially in VOI should intersect." );
    OPENGEODE_EXCEPTION(
        !SpatialDomainChecker< OwnerSegment2D >::intersects_domain(
            domain, seg_outside ),
        "[Test] Segment completely outside should not intersect." );
    OPENGEODE_EXCEPTION(
        SpatialDomainChecker< OwnerSegment2D >::intersects_domain(
            domain, seg_cross ),
        "[Test] Segment crossing the voi should intersect." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        test_spatial_domain_2D();

        geode::Logger::info( "TEST SpatialDomain SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}