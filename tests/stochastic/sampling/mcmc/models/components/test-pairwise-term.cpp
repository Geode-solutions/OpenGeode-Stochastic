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

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

geode::uuid init_object_set( geode::ObjectSets< geode::Point2D >& pattern )
{
    geode::Point2D p1{ { 0., 0. } };
    geode::Point2D p2{ { 1., 1. } };

    auto set_id = pattern.add_set( "default_name" );
    pattern.add_object( std::move( p1 ), set_id );
    pattern.add_object( std::move( p2 ), set_id );

    return set_id;
}

void test_normal_positive_pairwise( double gamma,
    const geode::ObjectSets< geode::Point2D >& pattern,
    const geode::uuid& set_id )
{
    auto interaction =
        std::make_unique< geode::EuclideanCutoffInteraction< geode::Point2D > >(
            2.1 );

    geode::PairwiseTerm< geode::Point2D > term(
        "strauss", gamma, { set_id }, std::move( interaction ) );
    auto neg_log_gamma = -std::log( gamma );

    // p1 and p2 interact → 1 pair
    auto total = term.total_log( pattern );

    OPENGEODE_EXCEPTION( total == neg_log_gamma * 1.,
        "[test pairwise] - total_log wrong value." );

    // Adding a third point close to p1 → expect new interactions
    geode::Point2D p3{ { 0.5, 0.5 } };

    geode::ObjectRef< geode::Point2D > p_ref{ p3, set_id };
    auto delta_add = term.delta_log_add( pattern, p_ref );
    // p3 interacts with p1 and p2 → 2 new pairs
    OPENGEODE_EXCEPTION( delta_add == neg_log_gamma * 2.,
        "[test pairwise] - delta_log_add wrong value." );
    geode::ObjectId obj_id{ 0, set_id };
    auto delta_remove = term.delta_log_remove( pattern, obj_id );
    // Removing p1 removes its interaction with p2 → 1 removed pair
    OPENGEODE_EXCEPTION( delta_remove == neg_log_gamma * -1.,
        "[test pairwise] - delta_log_remove wrong value." );

    auto delta_change = term.delta_log_change( pattern, obj_id, p_ref );
    // Replacing p1 with p3 changes interactions: p3 interacts with p2 → 1 pair
    // Old p1 interacted with p2 → 1 pair → no net change
    OPENGEODE_EXCEPTION(
        delta_change == 0., "[test pairwise] - delta_log_change wrong value." );
}

void test_zero_pairwise( double gamma,
    const geode::ObjectSets< geode::Point2D >& pattern,
    const geode::uuid& set_id )
{
    auto interaction =
        std::make_unique< geode::EuclideanCutoffInteraction< geode::Point2D > >(
            2.1 );

    geode::PairwiseTerm< geode::Point2D > term(
        "interaction", gamma, { set_id }, std::move( interaction ) );

    auto total = term.total_log( pattern );
    OPENGEODE_EXCEPTION(
        std::isinf( total ), "[test zero pairwise] - log_total wrong value." );

    geode::Point2D p3{ { 0.5, 0.5 } };
    geode::ObjectRef< geode::Point2D > p_ref{ p3, set_id };
    auto delta_add = term.delta_log_add( pattern, p_ref );
    OPENGEODE_EXCEPTION( std::isinf( delta_add ),
        "[test zero pairwise] - delta_log_add wrong value." );
    geode::ObjectId obj_id{ 0, set_id };
    auto delta_remove = term.delta_log_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION( delta_remove == 0.,
        "[test zero pairwise] - delta_log_remove wrong value." );

    auto delta_change = term.delta_log_change( pattern, obj_id, p_ref );
    OPENGEODE_EXCEPTION( delta_change == 0.,
        "[test zero pairwise] - delta_log_change wrong value." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        geode::ObjectSets< geode::Point2D > pattern;
        auto set_id = init_object_set( pattern );

        test_normal_positive_pairwise( 0.5, pattern, set_id );
        test_normal_positive_pairwise( geode::GLOBAL_EPSILON, pattern, set_id );
        test_normal_positive_pairwise( 2.0, pattern, set_id );

        test_zero_pairwise( 0., pattern, set_id );
        test_zero_pairwise( 0.9999 * geode::GLOBAL_EPSILON, pattern, set_id );
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}