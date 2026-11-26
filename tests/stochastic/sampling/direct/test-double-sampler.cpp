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
#include <numeric>
#include <optional>
#include <variant>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/distributions.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/distance.hpp>

const int NUMBER_OF_SAMPLES = 10000;

void test_double_sampler( geode::RandomEngine& engine )
{
    geode::UniformClosed< double > uniform_closed_double;
    double value = 1000.;
    value = geode::DoubleSampler::sample( engine, uniform_closed_double );
    geode::OpenGeodeStochasticStochasticException::test(
        value >= 0. && value <= 1., "[Uniform] -  value out of range." );

    geode::UniformClosedOpen< double > uniform_closedopen_double;
    value = 1000.;
    value = geode::DoubleSampler::sample( engine, uniform_closedopen_double );
    geode::OpenGeodeStochasticStochasticException::test(
        value >= 0. && value <= 1., "[Uniform] -  value out of range." );

    geode::TruncatedGaussian t_gaussian_double;
    t_gaussian_double.mean = 0.;
    t_gaussian_double.standard_deviation = 1.;
    t_gaussian_double.max_value = 1.;
    t_gaussian_double.min_value = -1.;
    value = 1000.;
    value = geode::DoubleSampler::sample( engine, t_gaussian_double );
    geode::OpenGeodeStochasticStochasticException::test(
        value >= -1. && value <= 1.,
        "[TruncatedGaussian] -  value out of range." );

    geode::Gaussian gaussian_double;
    t_gaussian_double.mean = 0.;
    t_gaussian_double.mean = 1.;

    value = 1000.;
    value = geode::DoubleSampler::sample( engine, t_gaussian_double );
    geode::OpenGeodeStochasticStochasticException::test(
        value != 1000., "[Gaussian] -  value did not changed." );
}

// Test function: create a uniform closed distribution from description
void test_create_uniform_closed( geode::RandomEngine& engine )
{
    geode::DoubleSampler::DistributionDescription desc;
    desc.name = "uniform_closed";
    desc.distribution_type =
        geode::UniformClosed< double >::distribution_type_static();
    desc.min_value = 2.;
    desc.max_value = 5.;

    auto dist = geode::DoubleSampler::create_distribution( desc );
    double value = geode::DoubleSampler::sample( engine, dist );
    geode::OpenGeodeStochasticStochasticException::test(
        value >= 2. && value <= 5., "[UniformClosed] - value out of bounds" );
}

// Test function: create a truncated Gaussian distribution from description
void test_create_truncated_gaussian( geode::RandomEngine& engine )
{
    geode::DoubleSampler::DistributionDescription desc;
    desc.name = "truncated_gaussian";
    desc.distribution_type =
        geode::TruncatedGaussian::distribution_type_static();
    desc.min_value = -2.;
    desc.max_value = 2.;
    desc.mean = 0.0;
    desc.standard_deviation = 1.0;

    auto dist = geode::DoubleSampler::create_distribution( desc );
    for( int i = 0; i < 100; ++i )
    {
        double value = geode::DoubleSampler::sample( engine, dist );
        geode::OpenGeodeStochasticStochasticException::test(
            value >= -2. && value <= 2.,
            "[TruncatedGaussian] - value out of bounds" );
    }
}

// An axial von Mises law is created from a description with mean and kappa,
// and samples an axis in [0, pi[
void test_create_axial_von_mises( geode::RandomEngine& engine )
{
    geode::DoubleSampler::DistributionDescription desc;
    desc.name = "axial_von_mises";
    desc.distribution_type = geode::AxialVonMises::distribution_type_static();
    desc.mean = -0.5;
    desc.kappa = 2.;

    const auto dist = geode::DoubleSampler::create_distribution( desc );
    const auto& law = std::get< geode::AxialVonMises >( dist );
    geode::OpenGeodeStochasticStochasticException::test(
        law.mean == -0.5 && law.concentration == 2.,
        "[AxialVonMises] - description parameters not transmitted" );
    for( int i = 0; i < 100; ++i )
    {
        const auto value = geode::DoubleSampler::sample( engine, dist );
        geode::OpenGeodeStochasticStochasticException::test(
            value >= 0. && value < M_PI,
            "[AxialVonMises] - value out of [0, pi[" );
    }

    desc.kappa.reset();
    bool rejected{ false };
    try
    {
        geode::DoubleSampler::create_distribution( desc );
    }
    catch( const geode::OpenGeodeException& )
    {
        rejected = true;
    }
    geode::OpenGeodeStochasticStochasticException::test(
        rejected, "[AxialVonMises] - description without kappa accepted" );
}

geode::DoubleSampler::Distribution angle_distribution(
    const geode::DistributionType& type,
    std::optional< double > min_value,
    std::optional< double > max_value,
    std::optional< double > mean,
    std::optional< double > standard_deviation )
{
    geode::DoubleSampler::DistributionDescription desc{ "angle" };
    desc.distribution_type = type;
    desc.min_value = min_value;
    desc.max_value = max_value;
    desc.mean = mean;
    desc.standard_deviation = standard_deviation;
    desc.kappa = 1.;
    return geode::DoubleSampler::create_angle_distribution_from_degrees( desc );
}

bool is_close( double value, double expected )
{
    return std::abs( value - expected ) < 1e-12;
}

// Degrees are converted to radians without normalizing the bounds nor the
// standard deviation.
void test_angle_from_degrees( geode::RandomEngine& engine )
{
    const auto full_circle =
        std::get< geode::UniformClosed< double > >( angle_distribution(
            geode::UniformClosed< double >::distribution_type_static(), 0.,
            360., std::nullopt, std::nullopt ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( full_circle.min_value, 0. )
            && is_close( full_circle.max_value, 2. * M_PI ),
        "[Angle] - 0-360 degrees should be [0, 2pi], got [",
        full_circle.min_value, ", ", full_circle.max_value, "]" );
    double max_sample{ 0. };
    for( const auto count : geode::Range{ 1000 } )
    {
        geode_unused( count );
        max_sample = std::max(
            max_sample, geode::DoubleSampler::sample( engine, full_circle ) );
    }
    geode::OpenGeodeStochasticStochasticException::test( max_sample > M_PI,
        "[Angle] - 0-360 degrees samples do not cover the circle" );

    // a closed-open range covers the circle once: 360 is never drawn
    const auto circle_once =
        std::get< geode::UniformClosedOpen< double > >( angle_distribution(
            geode::UniformClosedOpen< double >::distribution_type_static(), 0.,
            360., std::nullopt, std::nullopt ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( circle_once.min_value, 0. )
            && is_close( circle_once.max_value, 2. * M_PI ),
        "[Angle] - 0-360 degrees closed-open should be [0, 2pi[, got [",
        circle_once.min_value, ", ", circle_once.max_value, "[" );
    for( const auto count : geode::Range{ 1000 } )
    {
        geode_unused( count );
        const auto angle = geode::DoubleSampler::sample( engine, circle_once );
        geode::OpenGeodeStochasticStochasticException::test(
            angle >= 0. && angle < circle_once.max_value,
            "[Angle] - 0-360 degrees closed-open sample out of [0, 2pi[: ",
            angle );
    }

    const auto across_zero =
        std::get< geode::UniformClosed< double > >( angle_distribution(
            geode::UniformClosed< double >::distribution_type_static(), -90.,
            90., std::nullopt, std::nullopt ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( across_zero.min_value, -M_PI / 2. )
            && is_close( across_zero.max_value, M_PI / 2. ),
        "[Angle] - -90-90 degrees should be [-pi/2, pi/2], got [",
        across_zero.min_value, ", ", across_zero.max_value, "]" );
    const auto value = geode::DoubleSampler::sample( engine, across_zero );
    geode::OpenGeodeStochasticStochasticException::test(
        value >= -M_PI / 2. && value <= M_PI / 2.,
        "[Angle] - -90-90 degrees sample out of range" );

    const auto wide = std::get< geode::Gaussian >(
        angle_distribution( geode::Gaussian::distribution_type_static(),
            std::nullopt, std::nullopt, 0., 400. ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( wide.standard_deviation, 400. * M_PI / 180. ),
        "[Angle] - a standard deviation of 400 degrees became ",
        wide.standard_deviation );

    // the mean stays consistent with the bounds of a truncated law
    const auto truncated =
        std::get< geode::TruncatedGaussian >( angle_distribution(
            geode::TruncatedGaussian::distribution_type_static(), 350., 390.,
            370., 10. ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( truncated.mean, 370. * M_PI / 180. )
            && is_close( truncated.min_value.value(), 350. * M_PI / 180. )
            && is_close( truncated.max_value.value(), 390. * M_PI / 180. ),
        "[Angle] - truncated law around 360 degrees: mean ", truncated.mean );

    // the mean direction of a von Mises law is normalized in [0, 2pi[
    const auto direction = std::get< geode::VonMises >(
        angle_distribution( geode::VonMises::distribution_type_static(),
            std::nullopt, std::nullopt, 450., std::nullopt ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( direction.mean, M_PI / 2. ) && direction.concentration == 1.,
        "[Angle] - von Mises mean of 450 degrees became ", direction.mean );

    // the mean axis of an axial von Mises law is normalized in [0, pi[
    for( const auto& [degrees, expected] :
        { std::pair{ 190., 10. }, std::pair{ -30., 150. }, std::pair{ 0., 0. },
            std::pair{ 370., 10. } } )
    {
        const auto axis = std::get< geode::AxialVonMises >( angle_distribution(
            geode::AxialVonMises::distribution_type_static(), std::nullopt,
            std::nullopt, degrees, std::nullopt ) );
        geode::OpenGeodeStochasticStochasticException::test(
            is_close( axis.mean, expected * M_PI / 180. )
                && axis.concentration == 1.,
            "[Angle] - axial von Mises mean of ", degrees, " degrees became ",
            axis.mean );
    }

    // the log-normal parameters describe ln(x): converting x to radians
    // shifts the mean by ln(pi/180) and keeps the standard deviation
    const auto log_normal =
        std::get< geode::TruncatedLogNormal >( angle_distribution(
            geode::TruncatedLogNormal::distribution_type_static(), 10., 80.,
            std::log( 30. ), 0.5 ) );
    geode::OpenGeodeStochasticStochasticException::test(
        is_close( log_normal.mean, std::log( 30. * M_PI / 180. ) )
            && is_close( log_normal.standard_deviation, 0.5 )
            && is_close( log_normal.min_value.value(), 10. * M_PI / 180. )
            && is_close( log_normal.max_value.value(), 80. * M_PI / 180. ),
        "[Angle] - log-normal of median 30 degrees became ", log_normal.mean,
        ", ", log_normal.standard_deviation );
}

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::RandomEngine random_engine;

        test_double_sampler( random_engine );

        test_create_uniform_closed( random_engine );
        test_create_truncated_gaussian( random_engine );
        test_create_axial_von_mises( random_engine );
        test_angle_from_degrees( random_engine );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
