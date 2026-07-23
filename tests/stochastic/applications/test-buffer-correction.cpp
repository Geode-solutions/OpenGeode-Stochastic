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

#include <geode/stochastic/applications/poisson_process.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>

namespace
{
    constexpr auto burn_in_steps{ 2000 };
    constexpr auto realizations{ 2000 };
    constexpr auto nb_steps{ 200 };
    constexpr auto area_size{ 50. };
    constexpr auto lambda{ 0.05 };

    template < typename ObjectType >
    double run_poisson( geode::RandomEngine& engine,
        const geode::PoissonProcessBuilder< ObjectType >& poisson )
    {
        geode::SimulationRunner< ObjectType > runner{
            poisson.build_simulation_context()
        };

        auto& burn_state = runner.run( engine, burn_in_steps );
        geode_unused( burn_state );

        double mean_nb_object = 0.;
        for( const auto real_id : geode::Range{ realizations } )
        {
            geode_unused( real_id );
            auto& state = runner.run( engine, nb_steps );
            mean_nb_object += state.nb_objects();
        }

        return mean_nb_object / realizations;
    }

    void test_invariance( const std::vector< double >& means )
    {
        constexpr double tolerance = 0.001;
        for( const auto id : geode::Range{ means.size() - 1 } )
        {
            geode::OpenGeodeStochasticStochasticException::test(
                std::abs( means[id] - lambda ) < tolerance,
                "[Buffer] - Spatial buffer introduces a bias." );
        }
    }

    void test_point_buffer_invariance()
    {
        geode::Logger::info( "TEST - POINT PROCESS BUFFER INVARIANCE" );

        geode::RandomEngine engine;
        engine.set_seed( "@point-buffer@" );

        std::vector< double > means;

        // NOLINTBEGIN(*-magic-numbers)

        for( const double buffer : { 0.0, 10., 20.0, 50.0 } )
        {
            geode::SpatialDomainConfig< 2 > domain;
            domain.min_point = geode::Point2D{ { 0., 0. } };
            domain.max_point = geode::Point2D{ { area_size, area_size } };
            domain.buffer_size = buffer;

            const auto area = std::pow( area_size + ( 2.0 * buffer ), 2. );

            geode::PoissonProcessBuilder< geode::Point2D > poisson;
            poisson.set_domain( domain );

            auto& set_config =
                poisson.add_set( "points", lambda, std::nullopt );

            auto stats = run_poisson( engine, poisson );
            means.push_back( stats / area );
            geode::Logger::info( "buffer=", buffer,
                " mean density=", means.back(),
                " err= ", means.back() - lambda );
        }
        test_invariance( means );

        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_segment_buffer_invariance()
    {
        geode::Logger::info( "TEST - SEGMENT PROCESS BUFFER INVARIANCE" );

        geode::RandomEngine engine;
        engine.set_seed( "@segment-buffer@" );

        std::vector< double > means;

        // NOLINTBEGIN(*-magic-numbers)

        for( const double buffer : { 0.0, 10., 20.0, 50.0 } )
        {
            geode::SpatialDomainConfig< 2 > domain;
            domain.min_point = geode::Point2D{ { 0., 0. } };
            domain.max_point = geode::Point2D{ { area_size, area_size } };
            domain.buffer_size = buffer;

            const auto area = std::pow( area_size + ( 2.0 * buffer ), 2. );

            geode::PoissonProcessBuilder< geode::OwnerSegment2D > poisson;
            poisson.set_domain( domain );

            auto& set_config =
                poisson.add_set( "segments", lambda, std::nullopt );

            set_config.sampler.length.distribution_type =
                geode::UniformClosed< double >::distribution_type_static();
            set_config.sampler.length.min_value = 2;
            set_config.sampler.length.max_value = 5.;

            set_config.sampler.azimuth.distribution_type =
                geode::UniformClosed< double >::distribution_type_static();
            set_config.sampler.azimuth.min_value = 0;
            set_config.sampler.azimuth.max_value = 180.;

            auto stats = run_poisson( engine, poisson );
            means.push_back( stats / area );
            geode::Logger::info( "buffer=", buffer,
                " mean density=", means.back(),
                " err= ", means.back() - lambda );
        }
        test_invariance( means );

        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        test_point_buffer_invariance();
        test_segment_buffer_invariance();

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}