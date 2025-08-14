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

geode::Configuration< geode::Point2D > create_configuration()
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

void test_normal_positive_intensity(
    double lambda, const geode::Configuration< geode::Point2D > &pattern )
{
    geode::IntensityTerm< geode::Point2D > term( lambda );
    auto neg_log_lambda = -std::log( lambda );

    double total = term.log_total( pattern );
    OPENGEODE_EXCEPTION( total == neg_log_lambda * 2.,
        "[test intensity]- log_total wrong value." );

    geode::Point2D p3{ { 2., 2. } };
    geode::MarkedObject< geode::Point2D > mp3{ std::move( p3 ) };

    double delta_add = term.log_delta_add( pattern, mp3 );
    OPENGEODE_EXCEPTION( delta_add == neg_log_lambda * 1.,
        "[test intensity]- log_delta_add wrong value." );

    double delta_remove = term.log_delta_remove( pattern, 0 );

    OPENGEODE_EXCEPTION( delta_remove == neg_log_lambda * -1.,
        "[test intensity]- log_delta_remove wrong value." );

    double delta_change = term.log_delta_change( pattern, 0, mp3 );
    OPENGEODE_EXCEPTION(
        delta_change == 0., "[test intensity]- log_delta_change wrong value." );
}

void test_normal_zero_intensity(
    double lambda, const geode::Configuration< geode::Point2D > &pattern )
{
    geode::IntensityTerm< geode::Point2D > term( lambda );
    double total = term.log_total( pattern );

    OPENGEODE_EXCEPTION(
        std::isinf( total ), "[test zero intensity]- log_total wrong value." );

    geode::Point2D p3{ { 2., 2. } };
    geode::MarkedObject< geode::Point2D > mp3{ std::move( p3 ) };

    double delta_add = term.log_delta_add( pattern, mp3 );
    OPENGEODE_EXCEPTION( std::isinf( delta_add ),
        "[test zero intensity]- log_delta_add wrong value." );

    double delta_remove = term.log_delta_remove( pattern, 0 );
    OPENGEODE_EXCEPTION( delta_remove == 0.,
        "[test zero intensity]- log_delta_remove wrong value." );

    double delta_change = term.log_delta_change( pattern, 0, mp3 );
    OPENGEODE_EXCEPTION( delta_change == 0.,
        "[test zero intensity]- log_delta_change wrong value." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        auto pattern = create_configuration();

        test_normal_positive_intensity( 0.5, pattern );
        test_normal_positive_intensity( geode::GLOBAL_EPSILON, pattern );
        test_normal_positive_intensity( 100.0021165, pattern );

        test_normal_zero_intensity( 0., pattern );
        test_normal_zero_intensity( 0.9999 * geode::GLOBAL_EPSILON, pattern );
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }

    try
    {
        geode::StochasticLibrary::initialize();
        geode::IntensityTerm< geode::Point2D > term( -geode::GLOBAL_EPSILON );
        geode::Logger::info( "TEST FAILED" );
        return 1;
    }
    catch( geode::OpenGeodeException &expt )
    {
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
}