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

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_monitor.hpp>

#include <cmath>
#include <vector>

namespace
{
    void test_statistics_monitor_basic()
    {
        geode::Logger::info( "TEST - StatisticsMonitor basic functionality" );

        // --- Create monitor with 3 energy terms
        geode::StatisticsMonitor monitor( 3 );

        // --- Add 2 realizations
        monitor.add_realization( { 1.0, 2.0, 3.0 } );
        monitor.add_realization( { 2.0, 3.0, 4.0 } );

        // --- Check count
        geode::OpenGeodeStochasticStochasticException::test(
            monitor.statiscal_count() == 2, "Count mismatch" );

        const auto& means = monitor.means();
        const auto& variances = monitor.variances();

        // --- Check means
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( means[0] - 1.5 ) < 1e-12, "Mean[0] incorrect" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( means[1] - 2.5 ) < 1e-12, "Mean[1] incorrect" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( means[2] - 3.5 ) < 1e-12, "Mean[2] incorrect" );

        // --- Check variances
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( variances[0] - 0.5 ) < 1e-12, "Variance[0] incorrect" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( variances[1] - 0.5 ) < 1e-12, "Variance[1] incorrect" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( variances[2] - 0.5 ) < 1e-12, "Variance[2] incorrect" );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_statistics_monitor_basic();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}