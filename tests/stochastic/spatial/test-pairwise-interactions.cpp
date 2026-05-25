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

#include <geode/geometry/distance.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>
namespace
{
    void test_interaction()
    {
        geode::MinimalDistanceCutoff< geode::Point2D > interaction( 2.0 );

        geode::Point2D p1{ { 0., 0. } };
        geode::Point2D p2{ { 1., 1. } }; // distance = 1.44 < 2.0
        geode::uuid id{};

        double result =
            interaction.evaluate( geode::ObjectRef< geode::Point2D >{ p1, id },
                geode::ObjectRef< geode::Point2D >{ p2, id } );
        geode::OpenGeodeStochasticStochasticException::test( result == 1.0,
            "[MinimalDistanceCutoff] Failed for inside cutoff case." );
    }
    void test_no_interaction()
    {
        geode::MinimalDistanceCutoff< geode::Point2D > interaction( 2.0 );

        geode::Point2D p1{ { 0., 0. } };
        geode::Point2D p2{ { 3., 0. } }; // distance = 3.0 > 2.0
        geode::uuid id{};

        double result =
            interaction.evaluate( geode::ObjectRef< geode::Point2D >{ p1, id },
                geode::ObjectRef< geode::Point2D >{ p2, id } );
        geode::OpenGeodeStochasticStochasticException::test( result == 0.0,
            "[MinimalDistanceCutoff] Failed for outside cutoff case." );
    }
    void test_limit_interaction()
    {
        geode::MinimalDistanceCutoff< geode::Point2D > interaction( 2.0 );

        geode::Point2D p1{ { 0., 0. } };
        geode::Point2D p2{ { 2., 0. } }; // distance = 2.0 == cutoff
        geode::uuid id{};

        double result =
            interaction.evaluate( geode::ObjectRef< geode::Point2D >{ p1, id },
                geode::ObjectRef< geode::Point2D >{ p2, id } );
        geode::OpenGeodeStochasticStochasticException::test( result == 1.0,
            "[MinimalDistanceCutoff] Failed for exact cutoff case." );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();

        test_interaction();
        test_no_interaction();
        test_limit_interaction();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}