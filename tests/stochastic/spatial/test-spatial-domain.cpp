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
#include <geode/stochastic/spatial/single_object_features/single_object_feature.hpp>
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
    geode::OpenGeodeStochasticStochasticException::test(
        domain.n_volume() == 1.0, "Domain volume wrong." );
    geode::OpenGeodeStochasticStochasticException::test(
        domain.extended_n_volume() == 4,
        "Extended domain volume wrong." ); // (1+0.5*2)^2

    // Test points inside domain
    geode::Point2D inside{ { 0.5, 0.5 } };
    geode::Point2D boundary{ { 1., 1. } };
    geode::Point2D outside{ { 1.4, 1.4 } };

    geode::OpenGeodeStochasticStochasticException::test(
        domain.contains( inside ), "Point inside should be contained." );
    geode::OpenGeodeStochasticStochasticException::test(
        domain.contains( boundary ), "Boundary point should be contained." );
    geode::OpenGeodeStochasticStochasticException::test(
        !domain.contains( outside ), "Point outside should not be contained." );

    geode::OpenGeodeStochasticStochasticException::test(
        domain.extended_contains( outside ),
        "Point outside original but in buffer should be contained." );

    // Test SpatialDomainChecker for points: the closed domain
    geode::OpenGeodeStochasticStochasticException::test(
        SpatialDomainChecker< Point2D >::count_inside_domain( domain, inside )
            == 1.,
        "count_inside_domain failed for inside point" );
    geode::OpenGeodeStochasticStochasticException::test(
        SpatialDomainChecker< Point2D >::count_inside_domain( domain, boundary )
            == 1.,
        "count_inside_domain failed for boundary point" );
    geode::OpenGeodeStochasticStochasticException::test(
        SpatialDomainChecker< Point2D >::count_inside_domain( domain, outside )
            == 0.,
        "count_inside_domain failed for outside point" );
}

void check_feature(
    const geode::ObjectInDomainFeature< geode::OwnerSegment2D >& feature,
    const geode::SpatialDomain< 2 >& domain,
    const geode::OwnerSegment2D& segment,
    double expected,
    std::string_view label )
{
    const geode::OwnerSegment2D reversed{ segment.vertices()[1],
        segment.vertices()[0] };
    for( const auto* seg : { &segment, &reversed } )
    {
        const auto value = feature.evaluate_inside_domain( *seg, domain );
        geode::OpenGeodeStochasticStochasticException::test( value == expected,
            "[ObjectInDomainFeature] ", label, ": ", value, " instead of ",
            expected, " (reversed: ", seg == &reversed, ")" );
        geode::OpenGeodeStochasticStochasticException::test(
            feature.evaluate( *seg, domain ) == 1., "[ObjectInDomainFeature] ", label,
            ": a segment should count 1 in the extended domain" );
    }
}

// A segment is counted in the domain by its endpoints: each endpoint strictly
// inside counts 1/2, whichever is the anchor (point1).
void test_segment_in_domain_feature()
{
    const auto domain = init_domain();
    const geode::ObjectInDomainFeature< geode::OwnerSegment2D > feature;
    // NOLINTBEGIN(*-magic-numbers)
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { 0.1, 0.1 } }, geode::Point2D{ { 0.9, 0.9 } } },
        1., "both endpoints inside" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { 0.9, 0.9 } }, geode::Point2D{ { 1.4, 1.4 } } },
        0.5, "one endpoint inside" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { -0.4, 0.9 } }, geode::Point2D{ { 0.6, 0.6 } } },
        0.5, "one endpoint inside, the other in the buffer" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { -0.4, 0.5 } }, geode::Point2D{ { 1.4, 0.5 } } },
        0., "crossing the domain, no endpoint inside" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { 1.2, 1.2 } }, geode::Point2D{ { 1.4, 1.4 } } },
        0., "outside" );
    // A trace clipped by the window ends on its boundary: this cut point is
    // not an endpoint of the fracture, even with a rounding error.
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { 1., 0.5 } }, geode::Point2D{ { 0.5, 0.5 } } },
        0.5, "endpoint on the boundary" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{ geode::Point2D{ { 1. - 1e-12, 0.5 } },
            geode::Point2D{ { 0.5, 0.5 } } },
        0.5, "endpoint on the boundary up to rounding" );
    check_feature( feature, domain,
        geode::OwnerSegment2D{
            geode::Point2D{ { 0., 0. } }, geode::Point2D{ { 1., 1. } } },
        0., "diagonal clipped by the window" );
    // NOLINTEND(*-magic-numbers)
}

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        test_spatial_domain_2D();
        test_segment_in_domain_feature();

        geode::Logger::info( "TEST SpatialDomain SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}