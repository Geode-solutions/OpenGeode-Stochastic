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
#include <geode/stochastic/sampling/mcmc/energy_terms/density_term.hpp>

#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

geode::uuid init_object_set( geode::ObjectSets< geode::Point2D >& pattern )
{
    geode::Point2D p1{ { 0., 0. } };
    geode::Point2D p2{ { 1., 1. } };
    geode::Point2D p_buffer{ { 1.3, 0.1 } };

    auto set_id = pattern.add_set( "default_name" );
    pattern.add_object( std::move( p1 ), set_id, false );
    pattern.add_object( std::move( p2 ), set_id, false );
    pattern.add_object( std::move( p_buffer ), set_id, false ); // buffer last

    return set_id;
}

geode::SpatialDomain< 2 > init_domain()
{
    geode::BoundingBox2D box;
    box.add_point( geode::Point2D{ { -1.0, -1.0 } } );
    box.add_point( geode::Point2D{ { 1.0, 1.0 } } );
    return geode::SpatialDomain< 2 >{ box, 0.5 };
}

void run_density_test( double lambda,
    const geode::ObjectSets< geode::Point2D >& pattern,
    const geode::uuid& set_id,
    const geode::SpatialDomain< 2 >& domain )
{
    geode::DensityTerm< geode::Point2D > term(
        "density", lambda, { set_id }, domain );

    auto neg_log_lambda = -std::log( lambda );
    double expected_add =
        ( lambda > 0. ? neg_log_lambda
                      : std::numeric_limits< double >::infinity() );
    double expected_remove = ( lambda > 0. ? -neg_log_lambda : 0. );

    // --- Total log
    double expected_total =
        ( lambda > 0. ? neg_log_lambda * 2.
                      : std::numeric_limits< double >::infinity() );
    double total = term.total_log( pattern );
    OPENGEODE_EXCEPTION(
        total == expected_total, "[DensityTerm] total_log wrong" );

    // --- Delta add inside VOI
    geode::Point2D p_inside{ { 0.5, 0.5 } };
    geode::ObjectRef< geode::Point2D > ref_inside{ p_inside, set_id };
    double delta = term.delta_log_add( pattern, ref_inside );
    OPENGEODE_EXCEPTION(
        delta == expected_add, "[DensityTerm] delta_log_add inside VOI wrong" );

    // --- Delta add in buffer (outside VOI)
    geode::Point2D p_buffer{ { 1.3, 0.0 } };
    geode::ObjectRef< geode::Point2D > ref_buffer{ p_buffer, set_id };
    delta = term.delta_log_add( pattern, ref_buffer );
    OPENGEODE_EXCEPTION(
        delta == 0., "[DensityTerm] delta_log_add outside VOI wrong" );

    // --- Delta remove anchored object
    geode::ObjectId obj_id{ 0, false, set_id };
    delta = term.delta_log_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION(
        delta == expected_remove, "[DensityTerm] delta_log_remove wrong" );

    // --- Delta change anchored → buffer
    geode::ObjectRef< geode::Point2D > new_buffer{ p_buffer, set_id };
    delta = term.delta_log_change( pattern, obj_id, new_buffer );
    OPENGEODE_EXCEPTION( delta == expected_remove,
        "[DensityTerm] delta_log_change anchored→buffer wrong" );

    // --- Delta change anchored → anchored
    geode::Point2D p_anchored{ { 0.1, 0.1 } };
    geode::ObjectRef< geode::Point2D > new_anchored{ p_anchored, set_id };
    delta = term.delta_log_change( pattern, obj_id, new_anchored );
    OPENGEODE_EXCEPTION(
        delta == 0., "[DensityTerm] delta_log_change anchored→anchored wrong" );

    // --- Delta change buffer → anchored
    geode::ObjectId buffer_id{ 2, false, set_id };
    delta = term.delta_log_change( pattern, buffer_id, ref_inside );
    OPENGEODE_EXCEPTION( delta == expected_add,
        "[DensityTerm] delta_log_change buffer→anchored wrong" );
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

        // Test different lambda values including near-zero
        run_density_test( 0.5, pattern, set_id, domain );
        run_density_test( geode::GLOBAL_EPSILON, pattern, set_id, domain );
        run_density_test( 100.0021165, pattern, set_id, domain );
        run_density_test( 0., pattern, set_id, domain ); // zero lambda
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }

    try
    {
        geode::StochasticLibrary::initialize();
        geode::uuid set_id;
        auto domain = init_domain();

        geode::DensityTerm< geode::Point2D > term(
            "zero", -geode::GLOBAL_EPSILON, { set_id }, domain );
        geode::Logger::info( "TEST FAILED" );
        return 1;
    }
    catch( geode::OpenGeodeException& expt )
    {
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
}