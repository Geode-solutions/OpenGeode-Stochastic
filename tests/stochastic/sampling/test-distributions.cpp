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

#include <limits>

#include <geode/basic/logger.hpp>

#include <geode/stochastic/sampling/distributions.hpp>

// NOLINTBEGIN(*-magic-numbers)
namespace
{
    constexpr auto INF = std::numeric_limits< double >::infinity();

    template < typename Law >
    void check_valid( const Law& law, bool expected, std::string_view message )
    {
        geode::OpenGeodeStochasticStochasticException::test(
            law.is_valid() == expected, message );
    }

    void test_truncated_gaussian()
    {
        geode::TruncatedGaussian law;
        check_valid( law, true, "[Test] Default truncated gaussian is valid" );
        law.min_value = -1.;
        law.max_value = 2.;
        check_valid( law, true, "[Test] Bounded truncated gaussian is valid" );
        law.mean = 10.;
        law.standard_deviation = 2.;
        law.min_value = 5.;
        law.max_value = 15.;
        check_valid( law, true, "[Test] Shifted truncated gaussian is valid" );
        // range beyond 6 standard deviations of the mean
        law = {};
        law.min_value = 7.;
        check_valid( law, false, "[Test] Range above 6 sigma" );
        law.min_value = std::nullopt;
        law.max_value = -7.;
        check_valid( law, false, "[Test] Range below 6 sigma" );
        law.min_value = 5.5;
        law.max_value = 8.;
        check_valid( law, true, "[Test] Range overlapping the tail" );
        law = {};
        law.min_value = 2.;
        law.max_value = 1.;
        check_valid( law, false, "[Test] Inverted bounds" );
        law.max_value = 2.;
        check_valid( law, false, "[Test] Empty range" );
        // the 6 sigma limits are included
        law = {};
        law.min_value = 6.;
        check_valid( law, true, "[Test] Range starting at mean + 6 sigma" );
        law.min_value = std::nullopt;
        law.max_value = -6.;
        check_valid( law, true, "[Test] Range ending at mean - 6 sigma" );
        law.max_value = -6.1;
        check_valid( law, false, "[Test] Range ending below mean - 6 sigma" );
        law = {};
        law.standard_deviation = 0.;
        check_valid( law, false, "[Test] Null standard deviation" );
        law.standard_deviation = -1.;
        check_valid( law, false, "[Test] Negative standard deviation" );
        law.standard_deviation = INF;
        check_valid( law, false, "[Test] Infinite standard deviation" );
        law = {};
        law.mean = INF;
        check_valid( law, false, "[Test] Infinite mean" );
        law.mean = std::numeric_limits< double >::quiet_NaN();
        check_valid( law, false, "[Test] NaN mean" );
    }

    void test_truncated_lognormal()
    {
        geode::TruncatedLogNormal law;
        check_valid( law, true, "[Test] Default log normal is valid" );
        law.min_value = 0.5;
        law.max_value = 3.;
        check_valid( law, true, "[Test] Bounded log normal is valid" );
        law.max_value = std::nullopt;
        check_valid( law, true, "[Test] Lower bounded log normal is valid" );
        law.min_value = std::nullopt;
        law.max_value = 4.;
        check_valid( law, true, "[Test] Upper bounded log normal is valid" );
        law.min_value = 4.;
        check_valid( law, false, "[Test] Empty range" );
        law.min_value = 5.;
        check_valid( law, false, "[Test] Inverted bounds" );
        law.min_value = -1.;
        law.max_value = 4.;
        check_valid( law, false, "[Test] Negative lower bound" );
        law.min_value = 0.;
        check_valid( law, true, "[Test] Null lower bound" );
        // the lower bound is 0 when not given
        law.min_value = std::nullopt;
        law.max_value = 0.;
        check_valid( law, false, "[Test] Upper bound at the default lower" );
        law = {};
        law.standard_deviation = 0.;
        check_valid( law, false, "[Test] Null standard deviation" );
        law.standard_deviation = INF;
        check_valid( law, false, "[Test] Infinite standard deviation" );
        law = {};
        law.mean = INF;
        check_valid( law, false, "[Test] Infinite mean" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_truncated_gaussian();
        test_truncated_lognormal();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
// NOLINTEND(*-magic-numbers)
