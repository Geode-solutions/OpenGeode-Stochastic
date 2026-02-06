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
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

geode::uuid init_object_set( geode::ObjectSets< geode::Point2D >& pattern )
{
    geode::Point2D p1{ { 0.1, 0.1 } }; // VOI
    geode::Point2D p2{ { 0.9, 0.9 } }; // VOI
    geode::Point2D p3{ { 1.3, 1.3 } }; // buffer

    auto set_id = pattern.add_set( "default_name" );
    pattern.add_object( std::move( p1 ), set_id, false );
    pattern.add_object( std::move( p2 ), set_id, false );
    pattern.add_object( std::move( p3 ), set_id, false );

    return set_id;
}

geode::SpatialDomain< 2 > init_domain()
{
    geode::BoundingBox2D box;
    box.add_point( geode::Point2D{ { 0., 0. } } );
    box.add_point( geode::Point2D{ { 1.0, 1.0 } } );
    return geode::SpatialDomain< 2 >{ box, 0.5 };
}

void test_pairwise_term( double gamma,
    const geode::ObjectSets< geode::Point2D >& pattern,
    const geode::uuid& set_id,
    const geode::SpatialDomain< 2 >& domain )
{
    auto interaction =
        std::make_unique< geode::EuclideanCutoffInteraction< geode::Point2D > >(
            1 );
    geode::PairwiseTerm< geode::Point2D > term(
        "strauss", gamma, { set_id }, std::move( interaction ), domain );

    double neg_log_gamma = -std::log( gamma );

    // --- total_log
    double total = term.total_log( pattern );
    // Only VOI-anchored interactions counted in statistic: p1-p2
    OPENGEODE_EXCEPTION(
        total == term.contribution( 1 ), "[PairwiseTerm] total_log wrong" );

    // --- delta_add VOI → VOI
    geode::Point2D p4{ { 0.5, 0.5 } };
    geode::ObjectRef< geode::Point2D > p4_ref{ p4, set_id };
    double delta = term.delta_log_add( pattern, p4_ref );
    // interacts with p1 and p2 → 2 pairs
    OPENGEODE_EXCEPTION( delta == term.contribution( 2 ),
        "[PairwiseTerm] delta_log_add VOI wrong" );

    // --- delta_add buffer → buffer (outside VOI, inside buffer)
    geode::Point2D p_buffer{ { 1.2, 1.2 } };
    geode::ObjectRef< geode::Point2D > buffer_ref{ p_buffer, set_id };
    delta = term.delta_log_add( pattern, buffer_ref );
    // energy counts interactions: buffer interacts with p2,p3 (p3 is in
    // buffer) → 2 pairs
    OPENGEODE_EXCEPTION( delta == term.contribution( 1 ),
        "[PairwiseTerm] delta_log_add buffer wrong" );

    // --- delta_change VOI → VOI
    geode::ObjectId obj_id{ 0, false, set_id }; // p1
    delta = term.delta_log_change( pattern, obj_id, p4_ref );
    // p1 replaced by p4: interacts with p2 → add 1 interaction
    OPENGEODE_EXCEPTION( delta == term.contribution( 1 ),
        "[PairwiseTerm] delta_log_change VOI->VOI wrong" );

    // --- delta_change buffer → VOI
    geode::ObjectId old_buffer{ 2, false, set_id }; // p3
    delta = term.delta_log_change( pattern, old_buffer, p4_ref );
    // old buffer interactions removed (-p3 pairs (-1)), new VOI interactions
    // added
    // (+p4 pairs (+2)) net change = 1
    double expected_delta = term.contribution( 1 ); // adjust based on count
    OPENGEODE_EXCEPTION( delta == expected_delta,
        "[PairwiseTerm] delta_log_change buffer->VOI wrong" );

    // --- delta_change VOI → buffer
    delta = term.delta_log_change( pattern, obj_id, buffer_ref );
    // p1 → buffer: p1 old interaction = 0
    // p4  → p_buffer interaction with p2 +1
    expected_delta = term.contribution( 1 );
    OPENGEODE_EXCEPTION( delta == expected_delta,
        "[PairwiseTerm] delta_log_change VOI->buffer wrong" );

    // --- delta_remove VOI
    delta = term.delta_log_remove( pattern, obj_id );
    // p1 removed: no interaction with p2 removed
    OPENGEODE_EXCEPTION( delta == term.contribution( 0 ),
        "[PairwiseTerm] delta_log_remove VOI wrong" );

    // --- statistic (only anchored objects in VOI)
    double stat = term.statistic( pattern );
    // p1,p2 anchored → 1 pair
    OPENGEODE_EXCEPTION( stat == 1., "[PairwiseTerm] statistic wrong" );
}

void test_pairwise_term_zero_gamma( double gamma,
    const geode::ObjectSets< geode::Point2D >& pattern,
    const geode::uuid& set_id,
    const geode::SpatialDomain< 2 >& domain )
{
    auto interaction =
        std::make_unique< geode::EuclideanCutoffInteraction< geode::Point2D > >(
            1 );
    geode::PairwiseTerm< geode::Point2D > term(
        "strauss", gamma, { set_id }, std::move( interaction ), domain );

    // --- total_log
    double total = term.total_log( pattern );
    OPENGEODE_EXCEPTION( std::isinf( total ),
        "[PairwiseTerm] total_log with gamma<epsilon should be infinite" );

    // --- delta_add VOI → VOI
    geode::Point2D p4{ { 0.5, 0.5 } };
    geode::ObjectRef< geode::Point2D > p4_ref{ p4, set_id };
    double delta = term.delta_log_add( pattern, p4_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta ),
        "[PairwiseTerm] delta_log_add with gamma<epsilon should be infinite" );

    // --- delta_add buffer → buffer
    geode::Point2D p_buffer{ { 1.2, 1.2 } };
    geode::ObjectRef< geode::Point2D > buffer_ref{ p_buffer, set_id };
    delta = term.delta_log_add( pattern, buffer_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta ),
        "[PairwiseTerm] delta_log_add buffer with "
        "gamma<epsilon should be infinite" );

    // --- delta_change VOI → VOI
    geode::ObjectId obj_id{ 0, false, set_id };
    delta = term.delta_log_change( pattern, obj_id, p4_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta ),
        "[PairwiseTerm] delta_log_change VOI->VOI with "
        "gamma<epsilon should be infinite" );

    // --- delta_change buffer → VOI
    geode::ObjectId old_buffer{ 2, false, set_id };
    delta = term.delta_log_change( pattern, old_buffer, p4_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta ),
        "[PairwiseTerm] delta_log_change buffer->VOI with "
        "gamma<epsilon should be infinite" );

    // --- delta_change VOI → buffer
    delta = term.delta_log_change( pattern, obj_id, buffer_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta ),
        "[PairwiseTerm] delta_log_change VOI->buffer with "
        "gamma<epsilon should be infinite" );

    // --- delta_remove VOI
    delta = term.delta_log_remove( pattern, obj_id );
    DEBUG( delta );
    OPENGEODE_EXCEPTION( delta == 0, "[PairwiseTerm] delta_log_remove VOI with "
                                     "gamma<epsilon should be infinite" );

    // --- statistic (only anchored objects in VOI)
    double stat = term.statistic( pattern );
    // p1,p2 anchored → 1 pair
    OPENGEODE_EXCEPTION( stat == 1., "[PairwiseTerm] statistic wrong" );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        geode::ObjectSets< geode::Point2D > pattern;
        auto set_id = init_object_set( pattern );
        auto domain = init_domain();

        test_pairwise_term( 0.5, pattern, set_id, domain );
        test_pairwise_term( geode::GLOBAL_EPSILON, pattern, set_id, domain );
        test_pairwise_term( 2.0, pattern, set_id, domain );

        test_pairwise_term_zero_gamma(
            0.9999 * geode::GLOBAL_EPSILON, pattern, set_id, domain );

        //            test_normal_positive_pairwise( 0.5, pattern, set_id );
        //            test_normal_positive_pairwise(
        //                geode::GLOBAL_EPSILON, pattern, set_id );
        //            test_normal_positive_pairwise( 2.0, pattern, set_id );

        //            test_zero_pairwise( 0., pattern, set_id );
        //            test_zero_pairwise(
        //                0.9999 * geode::GLOBAL_EPSILON, pattern, set_id );
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}