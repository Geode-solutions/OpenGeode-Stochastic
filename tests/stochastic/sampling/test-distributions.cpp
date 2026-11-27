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
#include <limits>

#include <geode/basic/logger.hpp>

#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/distributions.hpp>

// NOLINTBEGIN(*-magic-numbers)
namespace
{
    constexpr auto INF = std::numeric_limits< double >::infinity();
    constexpr auto NaN = std::numeric_limits< double >::quiet_NaN();

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

    void test_truncated_powerlaw()
    {
        geode::TruncatedPowerLaw law;
        law.min_value = 1.;
        check_valid( law, true, "[Test] Lower bounded power law is valid" );
        law.max_value = 10.;
        check_valid( law, true, "[Test] Bounded power law is valid" );
        // x^-alpha with alpha <= 1 cannot be normalized up to infinity
        for( const auto alpha : { 1., 0.5 } )
        {
            law.alpha = alpha;
            law.max_value = 10.;
            check_valid( law, true, "[Test] Bounded power law, alpha <= 1" );
            law.max_value = std::nullopt;
            check_valid( law, false, "[Test] Unbounded power law, alpha <= 1" );
        }
        law = {};
        law.alpha = INF;
        check_valid( law, false, "[Test] Infinite exponent" );
        law = {};
        law.min_value = -1.;
        law.max_value = 2.;
        check_valid( law, false, "[Test] Negative lower bound" );
    }

    void test_invalid_description_is_rejected()
    {
        // create_distribution applies is_valid(): an invalid law fails when
        // the model is built, not at its first sample.
        geode::DoubleSampler::DistributionDescription description{ "law" };
        description.distribution_type =
            geode::TruncatedPowerLaw::distribution_type_static();
        description.alpha = 1.;
        description.min_value = 1.;
        bool rejected{ false };
        try
        {
            geode_unused(
                geode::DoubleSampler::create_distribution( description ) );
        }
        catch( const geode::OpenGeodeException& exception )
        {
            rejected =
                exception.type() == geode::OpenGeodeException::TYPE::data;
        }
        geode::OpenGeodeStochasticStochasticException::test(
            rejected, "[Test] Invalid description accepted" );
        description.max_value = 10.;
        geode_unused(
            geode::DoubleSampler::create_distribution( description ) );
    }

    void test_uniform()
    {
        geode::UniformClosed< double > closed;
        check_valid( closed, true, "[Test] Default uniform is valid" );
        closed.max_value = 0.;
        check_valid( closed, true, "[Test] Equal bounds are valid" );
        closed.max_value = INF;
        check_valid( closed, false, "[Test] Infinite upper bound" );
        closed.min_value = -INF;
        closed.max_value = 1.;
        check_valid( closed, false, "[Test] Infinite lower bound" );
        closed.min_value = 0.;
        closed.max_value = NaN;
        check_valid( closed, false, "[Test] NaN upper bound" );
        geode::UniformClosedOpen< double > closed_open;
        check_valid( closed_open, true, "[Test] Default uniform is valid" );
        closed_open.max_value = INF;
        check_valid( closed_open, false, "[Test] Infinite upper bound" );
        closed_open.max_value = NaN;
        check_valid( closed_open, false, "[Test] NaN upper bound" );
        closed_open = {};
        geode::OpenGeodeStochasticStochasticException::test(
            closed_open.string() == "UniformClosedOpen[0,1)",
            "[Test] Wrong closed open interval display: ",
            closed_open.string() );
        closed = {};
        geode::OpenGeodeStochasticStochasticException::test(
            closed.string() == "UniformClosed[0,1]",
            "[Test] Wrong closed interval display: ", closed.string() );
        geode::UniformClosed< geode::index_t > indices;
        indices.max_value = std::numeric_limits< geode::index_t >::max();
        check_valid( indices, true, "[Test] Integer uniform is valid" );
    }

    void test_nan_bounds()
    {
        // a comparison with NaN is always false: NaN bounds must be rejected
        geode::TruncatedGaussian gaussian;
        gaussian.min_value = NaN;
        check_valid( gaussian, false, "[Test] Gaussian NaN lower bound" );
        gaussian = {};
        gaussian.max_value = NaN;
        check_valid( gaussian, false, "[Test] Gaussian NaN upper bound" );
        geode::TruncatedLogNormal lognormal;
        lognormal.min_value = NaN;
        check_valid( lognormal, false, "[Test] Log normal NaN lower bound" );
        lognormal = {};
        lognormal.max_value = NaN;
        check_valid( lognormal, false, "[Test] Log normal NaN upper bound" );
        geode::TruncatedPowerLaw power_law;
        power_law.min_value = NaN;
        power_law.max_value = 2.;
        check_valid( power_law, false, "[Test] Power law NaN lower bound" );
        power_law.min_value = 1.;
        power_law.max_value = NaN;
        check_valid( power_law, false, "[Test] Power law NaN upper bound" );
    }

    void test_von_mises()
    {
        geode::VonMises law;
        check_valid( law, true, "[Test] Default von Mises is valid" );
        law.concentration = 0.;
        check_valid( law, true, "[Test] Null concentration is valid" );
        law.concentration = INF;
        check_valid( law, false, "[Test] Infinite concentration" );
        law.concentration = -1.;
        check_valid( law, false, "[Test] Negative concentration" );
        law.concentration = 1.;
        law.mean = NaN;
        check_valid( law, false, "[Test] NaN mean" );
    }

    void test_axial_von_mises()
    {
        geode::AxialVonMises law;
        check_valid( law, true, "[Test] Default axial von Mises is valid" );
        law.concentration = 0.;
        check_valid( law, true, "[Test] Null axial concentration is valid" );
        law.concentration = INF;
        check_valid( law, false, "[Test] Infinite axial concentration" );
        law.concentration = -1.;
        check_valid( law, false, "[Test] Negative axial concentration" );
        law.concentration = 1.;
        law.mean = NaN;
        check_valid( law, false, "[Test] NaN axial mean" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_truncated_gaussian();
        test_truncated_lognormal();
        test_truncated_powerlaw();
        test_uniform();
        test_nan_bounds();
        test_von_mises();
        test_axial_von_mises();
        test_invalid_description_is_rejected();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
// NOLINTEND(*-magic-numbers)
