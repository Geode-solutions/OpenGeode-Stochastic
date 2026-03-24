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
#include <geode/stochastic/sampling/mcmc/energy_terms/intensity_term.hpp>

#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

geode::uuid init_segment_set(
    geode::ObjectSets< geode::OwnerSegment2D >& pattern )
{
    using Segment = geode::OwnerSegment2D;

    // Fully inside domain
    Segment s1{ geode::Point2D{ { -0.5, 0.0 } },
        geode::Point2D{ { 0.5, 0.0 } } }; // length = 1.0

    // Partially inside domain
    Segment s2{ geode::Point2D{ { -2.0, 0.0 } },
        geode::Point2D{ { 0.0, 0.0 } } }; // clipped length = 1.0

    // Fully outside domain (buffer)
    Segment s_buffer{ geode::Point2D{ { 2.0, 2.0 } },
        geode::Point2D{ { 3.0, 3.0 } } }; // clipped = 0.0

    auto set_id = pattern.add_set( "segments" );
    pattern.add_object( std::move( s1 ), set_id, false );
    pattern.add_object( std::move( s2 ), set_id, false );
    pattern.add_object( std::move( s_buffer ), set_id, false );

    return set_id;
}

geode::SpatialDomain< 2 > init_domain()
{
    geode::BoundingBox2D box;
    box.add_point( geode::Point2D{ { -1.0, -1.0 } } );
    box.add_point( geode::Point2D{ { 1.0, 1.0 } } );
    return geode::SpatialDomain< 2 >{ box, 0.5 };
}

void run_intensity_test( double lambda,
    double characteristic_length,
    const geode::ObjectSets< geode::OwnerSegment2D >& pattern,
    const geode::uuid& set_id,
    const geode::SpatialDomain< 2 >& domain )
{
    geode::IntensityTerm term(
        "intensity", lambda, { set_id }, characteristic_length, domain );

    const double neg_log_lambda =
        ( lambda > 0. ? -std::log( lambda )
                      : std::numeric_limits< double >::infinity() );

    // Total clipped length inside domain:
    // s1: 1.0
    // s2: 1.0
    // s_buffer: 0.0
    const double total_length = 2.0;
    const double scaled_total = total_length / characteristic_length;

    const double expected_total =
        ( lambda > 0. ? neg_log_lambda * scaled_total
                      : std::numeric_limits< double >::infinity() );

    // --- Total log
    double total = term.total_log( pattern );
    OPENGEODE_EXCEPTION(
        total == expected_total, "[IntensityTerm] total_log wrong" );

    // --- Delta add (segment fully inside)
    geode::OwnerSegment2D s_inside{ geode::Point2D{ { 0.0, -0.5 } },
        geode::Point2D{ { 0.0, 0.5 } } }; // length 1
    geode::ObjectRef< geode::OwnerSegment2D > ref_inside{ s_inside, set_id };

    double expected_add =
        ( lambda > 0. ? neg_log_lambda * ( 1.0 / characteristic_length )
                      : std::numeric_limits< double >::infinity() );

    double delta = term.delta_log_add( pattern, ref_inside );
    OPENGEODE_EXCEPTION(
        delta == expected_add, "[IntensityTerm] delta_log_add inside wrong" );

    // --- Delta add outside domain
    geode::OwnerSegment2D s_outside{ geode::Point2D{ { 2.0, 2.0 } },
        geode::Point2D{ { 3.0, 3.0 } } };
    geode::ObjectRef< geode::OwnerSegment2D > ref_out{ s_outside, set_id };

    delta = term.delta_log_add( pattern, ref_out );
    OPENGEODE_EXCEPTION(
        delta == 0.0, "[IntensityTerm] delta_log_add outside wrong" );

    // --- Delta remove (first segment)
    geode::ObjectId obj_id{ 0, false, set_id };
    double expected_remove = ( lambda > 0. ? -expected_add : 0.0 );

    delta = term.delta_log_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION(
        delta == expected_remove, "[IntensityTerm] delta_log_remove wrong" );

    // --- Delta change: inside → outside
    delta = term.delta_log_change( pattern, obj_id, ref_out );
    OPENGEODE_EXCEPTION( delta == expected_remove,
        "[IntensityTerm] delta_log_change inside→outside wrong" );

    // --- Delta change: outside → inside
    geode::ObjectId buffer_id{ 2, false, set_id };
    delta = term.delta_log_change( pattern, buffer_id, ref_inside );
    OPENGEODE_EXCEPTION( delta == expected_add,
        "[IntensityTerm] delta_log_change outside→inside wrong" );

    // --- Delta change: inside → inside (same length)
    geode::OwnerSegment2D s_same{ geode::Point2D{ { -0.2, 0.0 } },
        geode::Point2D{ { 0.8, 0.0 } } }; // length still 1
    geode::ObjectRef< geode::OwnerSegment2D > ref_same{ s_same, set_id };

    delta = term.delta_log_change( pattern, obj_id, ref_same );
    OPENGEODE_EXCEPTION(
        delta == 0.0, "[IntensityTerm] delta_log_change inside→inside wrong" );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        geode::ObjectSets< geode::OwnerSegment2D > pattern;
        auto set_id = init_segment_set( pattern );
        auto domain = init_domain();

        const double L0 = 1.0;

        run_intensity_test( 0.5, L0, pattern, set_id, domain );
        run_intensity_test(
            geode::GLOBAL_EPSILON, L0, pattern, set_id, domain );
        run_intensity_test( 50.0, L0, pattern, set_id, domain );
        run_intensity_test( 0.0, L0, pattern, set_id, domain );
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}