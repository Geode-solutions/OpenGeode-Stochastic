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
#include <geode/stochastic/sampling/mcmc/models/components/intensity_term.hpp>

#include <geode/geometry/point.hpp>
#include <geode/stochastic/configuration/configuration.hpp>

geode::Configuration< geode::Point2D > create_configuration(
    const geode::uuid &group_id )
{
    geode::Point2D p1{ { 0., 0. } };
    geode::Point2D p2{ { 1., 1. } };

    geode::Configuration< geode::Point2D > pattern;
    pattern.add_group( group_id );
    pattern.add_object( std::move( p1 ), group_id );
    pattern.add_object( std::move( p2 ), group_id );

    return pattern;
}

void test_normal_positive_intensity( double lambda,
    const geode::Configuration< geode::Point2D > &pattern,
    const geode::uuid &group_id )
{
    geode::IntensityTerm< geode::Point2D > term( lambda, group_id );
    auto neg_log_lambda = -std::log( lambda );

    double total = term.total_log( pattern );
    OPENGEODE_EXCEPTION( total == neg_log_lambda * 2.,
        "[test intensity]- total_log wrong value." );

    geode::Point2D p3{ { 2., 2. } };

    double delta_add = term.delta_log_add( pattern, p3, group_id );
    OPENGEODE_EXCEPTION( delta_add == neg_log_lambda * 1.,
        "[test intensity]- delta_log_add wrong value." );

    geode::ObjectId obj_id{ 0, group_id };
    double delta_remove = term.delta_log_remove( pattern, obj_id );

    OPENGEODE_EXCEPTION( delta_remove == neg_log_lambda * -1.,
        "[test intensity]- delta_log_remove wrong value." );

    double delta_change = term.delta_log_change( pattern, obj_id, p3 );
    OPENGEODE_EXCEPTION(
        delta_change == 0., "[test intensity]- delta_log_change wrong value." );
}

void test_normal_zero_intensity( double lambda,
    const geode::Configuration< geode::Point2D > &pattern,
    const geode::uuid &group_id )
{
    geode::IntensityTerm< geode::Point2D > term( lambda, group_id );
    double total = term.total_log( pattern );

    OPENGEODE_EXCEPTION(
        std::isinf( total ), "[test zero intensity]- total_log wrong value." );

    geode::Point2D p3{ { 2., 2. } };

    double delta_add = term.delta_log_add( pattern, p3, group_id );
    OPENGEODE_EXCEPTION( std::isinf( delta_add ),
        "[test zero intensity]- delta_log_add wrong value." );

    geode::ObjectId obj_id{ 0, group_id };

    double delta_remove = term.delta_log_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION( delta_remove == 0.,
        "[test zero intensity]- delta_log_remove wrong value." );

    double delta_change = term.delta_log_change( pattern, obj_id, p3 );
    OPENGEODE_EXCEPTION( delta_change == 0.,
        "[test zero intensity]- delta_log_change wrong value." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::uuid group_id;

        auto pattern = create_configuration( group_id );

        test_normal_positive_intensity( 0.5, pattern, group_id );
        test_normal_positive_intensity(
            geode::GLOBAL_EPSILON, pattern, group_id );
        test_normal_positive_intensity( 100.0021165, pattern, group_id );

        test_normal_zero_intensity( 0., pattern, group_id );
        test_normal_zero_intensity(
            0.9999 * geode::GLOBAL_EPSILON, pattern, group_id );
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }

    try
    {
        geode::StochasticLibrary::initialize();
        geode::uuid group_id;
        geode::IntensityTerm< geode::Point2D > term(
            -geode::GLOBAL_EPSILON, group_id );
        geode::Logger::info( "TEST FAILED" );
        return 1;
    }
    catch( geode::OpenGeodeException &expt )
    {
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
}