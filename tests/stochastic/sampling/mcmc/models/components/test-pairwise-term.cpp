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

#include <geode/geometry/point.hpp>
#include <geode/stochastic/configuration/configuration.hpp>

geode::Configuration< geode::Point2D > create_pairwise_configuration()
{
    geode::Point2D p1{ { 0., 0. } };
    geode::MarkedObject< geode::Point2D > mp1{ std::move( p1 ) };
    geode::Point2D p2{ { 1., 1. } };
    geode::MarkedObject< geode::Point2D > mp2{ std::move( p2 ) };

    geode::Configuration< geode::Point2D > pattern;
    pattern.add_object( std::move( mp1 ) );
    pattern.add_object( std::move( mp2 ) );

    return pattern;
}

void test_normal_positive_pairwise(
    double gamma, const geode::Configuration< geode::Point2D >& pattern )
{
    auto interaction_fn = []( const geode::MarkedObject< geode::Point2D >& a,
                              const geode::MarkedObject< geode::Point2D >& b ) {
        // Interaction if distance < sqrt(2.1) for example
        double dx = a.geometry().value( 0 ) - b.geometry().value( 0 );
        double dy = a.geometry().value( 1 ) - b.geometry().value( 1 );
        double dist_sq = dx * dx + dy * dy;
        return dist_sq < 2.1;
    };

    geode::PairwiseTerm< geode::Point2D, decltype( interaction_fn ) > term(
        gamma, interaction_fn );
    auto neg_log_gamma = -std::log( gamma );

    // p1 and p2 interact → 1 pair
    double total = term.log_total( pattern );
    OPENGEODE_EXCEPTION( total == neg_log_gamma * 1.,
        "[test pairwise] - log_total wrong value." );

    // Adding a third point close to p1 → expect new interactions
    geode::Point2D p3{ { 0.5, 0.5 } };
    geode::MarkedObject< geode::Point2D > mp3{ std::move( p3 ) };

    double delta_add = term.log_delta_add( pattern, mp3 );
    // p3 interacts with p1 and p2 → 2 new pairs
    OPENGEODE_EXCEPTION( delta_add == neg_log_gamma * 2.,
        "[test pairwise] - log_delta_add wrong value." );

    double delta_remove = term.log_delta_remove( pattern, 0 );
    // Removing p1 removes its interaction with p2 → 1 removed pair
    OPENGEODE_EXCEPTION( delta_remove == neg_log_gamma * -1.,
        "[test pairwise] - log_delta_remove wrong value." );

    double delta_change = term.log_delta_change( pattern, 0, mp3 );
    // Replacing p1 with p3 changes interactions: p3 interacts with p2 → 1 pair
    // Old p1 interacted with p2 → 1 pair → no net change
    OPENGEODE_EXCEPTION(
        delta_change == 0., "[test pairwise] - log_delta_change wrong value." );
}

void test_zero_pairwise(
    double gamma, const geode::Configuration< geode::Point2D >& pattern )
{
    auto interaction_fn = []( const geode::MarkedObject< geode::Point2D >& a,
                              const geode::MarkedObject< geode::Point2D >& b ) {
        double dx = a.geometry().value( 0 ) - b.geometry().value( 0 );
        double dy = a.geometry().value( 1 ) - b.geometry().value( 1 );
        double dist_sq = dx * dx + dy * dy;
        return dist_sq < 2.1;
    };

    geode::PairwiseTerm< geode::Point2D, decltype( interaction_fn ) > term(
        gamma, interaction_fn );

    double total = term.log_total( pattern );
    OPENGEODE_EXCEPTION(
        std::isinf( total ), "[test zero pairwise] - log_total wrong value." );

    geode::Point2D p3{ { 0.5, 0.5 } };
    geode::MarkedObject< geode::Point2D > mp3{ std::move( p3 ) };

    double delta_add = term.log_delta_add( pattern, mp3 );
    OPENGEODE_EXCEPTION( std::isinf( delta_add ),
        "[test zero pairwise] - log_delta_add wrong value." );

    double delta_remove = term.log_delta_remove( pattern, 0 );
    OPENGEODE_EXCEPTION( delta_remove == 0.,
        "[test zero pairwise] - log_delta_remove wrong value." );

    double delta_change = term.log_delta_change( pattern, 0, mp3 );
    OPENGEODE_EXCEPTION( delta_change == 0.,
        "[test zero pairwise] - log_delta_change wrong value." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        auto pattern = create_pairwise_configuration();

        test_normal_positive_pairwise( 0.5, pattern );
        test_normal_positive_pairwise( geode::GLOBAL_EPSILON, pattern );
        test_normal_positive_pairwise( 2.0, pattern );

        test_zero_pairwise( 0., pattern );
        test_zero_pairwise( 0.9999 * geode::GLOBAL_EPSILON, pattern );
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}