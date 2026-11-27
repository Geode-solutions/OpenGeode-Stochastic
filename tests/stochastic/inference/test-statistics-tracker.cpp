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
#include <cmath>
#include <vector>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/applications/strauss_process.hpp>
#include <geode/stochastic/inference/statistics_tracker.hpp>

namespace
{
    // Model with two terms: the density of the set and the interaction
    geode::SimulationContext< geode::Point2D > two_term_context()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { 10., 10. } };
        domain.buffer_size = 1.;
        geode::StraussProcessBuilder< geode::Point2D > strauss;
        strauss.set_domain( domain );
        auto& set = strauss.add_set( "points", 0.5 );
        geode_unused( set );
        strauss.add_interaction( { "points" }, 0.5, 1., std::nullopt );
        // NOLINTEND(*-magic-numbers)
        return strauss.build_simulation_context();
    }

    bool is_close( double value, double expected )
    {
        constexpr double tolerance{ 1e-12 };
        return std::fabs( value - expected ) < tolerance;
    }

    void test_means_and_variances()
    {
        geode::Logger::info( "TEST - STATISTICS TRACKER MEANS AND VARIANCES" );
        const auto context = two_term_context();
        const auto& model = *context.model;
        geode::StatisticsTracker< geode::Point2D > tracker{ model };
        geode::OpenGeodeStochasticStochasticException::test(
            tracker.nb_realizations() == 0 && tracker.variances()[0] == 0.,
            "[StatisticsTracker] An empty tracker should have no realization "
            "and a zero variance." );

        // NOLINTBEGIN(*-magic-numbers)
        tracker.add_realization( { 1., 2. } );
        geode::OpenGeodeStochasticStochasticException::test(
            tracker.variances()[0] == 0.,
            "[StatisticsTracker] The variance of one realization should be "
            "zero." );
        tracker.add_realization( { 2., 6. } );
        tracker.add_realization( { 3., 7. } );

        geode::OpenGeodeStochasticStochasticException::test(
            tracker.nb_realizations() == 3,
            "[StatisticsTracker] Wrong number of realizations." );
        const auto& means = tracker.means();
        const auto variances = tracker.variances();
        geode::OpenGeodeStochasticStochasticException::test(
            is_close( means[0], 2. ) && is_close( means[1], 5. ),
            "[StatisticsTracker] Wrong means." );
        // unbiased variance (divided by n - 1)
        geode::OpenGeodeStochasticStochasticException::test(
            is_close( variances[0], 1. ) && is_close( variances[1], 7. ),
            "[StatisticsTracker] Wrong variances." );

        // access by term, in the order of the model terms
        const auto names = tracker.term_names();
        geode::OpenGeodeStochasticStochasticException::test(
            names == model.term_names() && names.size() == 2,
            "[StatisticsTracker] The term names should be the model ones." );
        const auto interaction = model.terms().get_term_uuid( names[1] );
        geode::OpenGeodeStochasticStochasticException::test(
            is_close( tracker.mean( interaction ), 5. )
                && is_close( tracker.variance( interaction ), 7. ),
            "[StatisticsTracker] Wrong mean or variance of a term." );
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    // Welford's update: the variance of values with a large mean is exact,
    // where the sum of the squares would lose all the significant digits.
    void test_large_mean_variance()
    {
        geode::Logger::info( "TEST - STATISTICS TRACKER WITH A LARGE MEAN" );
        const auto context = two_term_context();
        geode::StatisticsTracker< geode::Point2D > tracker{ *context.model };
        // NOLINTBEGIN(*-magic-numbers)
        constexpr double offset{ 1e9 };
        for( const double value : { 1., 2., 3. } )
        {
            tracker.add_realization( { offset + value, -offset - value } );
        }
        const auto variances = tracker.variances();
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( variances[0] - 1. ) < 1e-6
                && std::fabs( variances[1] - 1. ) < 1e-6,
            "[StatisticsTracker] The variance should not depend on the mean: ",
            variances[0], " ", variances[1] );
        geode::OpenGeodeStochasticStochasticException::test(
            tracker.means()[0] == offset + 2.
                && tracker.means()[1] == -offset - 2.,
            "[StatisticsTracker] Wrong means with a large offset." );
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    // A realization with more or fewer values than the model terms is
    // refused, and leaves the tracker unchanged.
    void test_wrong_number_of_values()
    {
        geode::Logger::info(
            "TEST - STATISTICS TRACKER WRONG NUMBER OF VALUES" );
        const auto context = two_term_context();
        geode::StatisticsTracker< geode::Point2D > tracker{ *context.model };
        // NOLINTBEGIN(*-magic-numbers)
        tracker.add_realization( { 1., 2. } );
        for( const auto& values :
            std::vector< std::vector< double > >{ { 1. }, { 1., 2., 3. }, {} } )
        {
            bool rejected{ false };
            try
            {
                tracker.add_realization( values );
            }
            catch( const geode::OpenGeodeException& )
            {
                rejected = true;
            }
            geode::OpenGeodeStochasticStochasticException::test( rejected,
                "[StatisticsTracker] A realization with ", values.size(),
                " values should be refused for 2 terms." );
        }
        geode::OpenGeodeStochasticStochasticException::test(
            tracker.nb_realizations() == 1 && tracker.means()[0] == 1.
                && tracker.means()[1] == 2.,
            "[StatisticsTracker] A refused realization should not change the "
            "tracker." );
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_means_and_variances();
        test_large_mean_variance();
        test_wrong_number_of_values();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
