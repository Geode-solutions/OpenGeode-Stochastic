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

#include <geode/geometry/distance.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/configuration/marked_object.hpp>
#include <geode/stochastic/configuration/pairwise_interactions.hpp>

namespace geode
{
    void test_distance_cutoff_interaction()
    {
        // Cutoff distance = 2.0
        DistanceCutoffInteraction interaction( 2.0 );

        // Case 1: Inside cutoff
        {
            geode::Point2D p1{ { 0., 0. } };
            geode::MarkedObject< geode::Point2D > mp1{ std::move( p1 ) };
            geode::Point2D p2{ { 1., 1. } }; // distance = 1.44 < 2.0
            geode::MarkedObject< geode::Point2D > mp2{ std::move( p2 ) };

            double result = interaction( mp1, mp2 );
            OPENGEODE_EXCEPTION( result == 1.0,
                "[DistanceCutoffInteraction] Failed for inside cutoff case." );
        }

        // Case 2: Outside cutoff
        {
            geode::Point2D p1{ { 0., 0. } };
            geode::MarkedObject< geode::Point2D > mp1{ std::move( p1 ) };
            geode::Point2D p2{ { 3., 0. } }; // distance = 3.0 > 2.0
            geode::MarkedObject< geode::Point2D > mp2{ std::move( p2 ) };

            double result = interaction( mp1, mp2 );
            OPENGEODE_EXCEPTION( result == 0.0,
                "[DistanceCutoffInteraction] Failed for outside cutoff case." );
        }

        // Case 3: Exactly at cutoff
        {
            geode::Point2D p1{ { 0., 0. } };
            geode::MarkedObject< geode::Point2D > mp1{ std::move( p1 ) };
            geode::Point2D p2{ { 2., 0. } }; // distance = 2.0 == cutoff
            geode::MarkedObject< geode::Point2D > mp2{ std::move( p2 ) };

            double result = interaction( mp1, mp2 );
            OPENGEODE_EXCEPTION( result == 1.0,
                "[DistanceCutoffInteraction] Failed for exact cutoff case." );
        }
    }
} // namespace geode

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        geode::test_distance_cutoff_interaction();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}