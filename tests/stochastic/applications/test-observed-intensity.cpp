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

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/applications/observed_intensity.hpp>
#include <geode/stochastic/applications/poisson_process.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace
{
    // NOLINTBEGIN(*-magic-numbers)
    geode::SpatialDomainConfig< 2 > square_domain( double size, double buffer )
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { size, size } };
        domain.buffer_size = buffer;
        return domain;
    }

    void set_sampler( geode::ObjectSetDefinition< geode::Fracture >& fset,
        double min_length,
        double max_length,
        double min_azimuth, // radians
        double max_azimuth )
    {
        fset.sampler.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.length.min_value = min_length;
        fset.sampler.length.max_value = max_length;
        fset.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.azimuth.min_value = min_azimuth;
        fset.sampler.azimuth.max_value = max_azimuth;
    }

    void check_area(
        double area, double expected, double tolerance, std::string_view label )
    {
        geode::Logger::info( "Touching area (", label, "): ", area,
            " (expected: ", expected, ")" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( area - expected ) < tolerance, "[ObservedIntensity] ",
            label, ": touching area ", area, " instead of ", expected );
    }

    // Area of the anchors for which a segment touches a w x h domain:
    // w * h + E[L * (w |cos theta| + h |sin theta|)]
    void test_touching_domain_area()
    {
        geode::Logger::info( "TEST - TOUCHING DOMAIN AREA" );
        const auto domain_cfg = square_domain( 10., 5. );
        const auto domain = geode::build_spatial_domain( domain_cfg );

        geode::ObjectSetDefinition< geode::Fracture > fset;
        set_sampler( fset, 4., 4., M_PI / 6., M_PI / 6. );
        const auto azimuth = M_PI / 6.;
        check_area( geode::touching_domain_area( *domain, fset.sampler ),
            100.
                + 4.
                      * ( 10. * std::abs( std::cos( azimuth ) )
                          + 10. * std::abs( std::sin( azimuth ) ) ),
            1e-9, "length 4, azimuth pi / 6" );

        set_sampler( fset, 2., 4., 0., M_PI );
        // E[L] = 3, E|cos| = E|sin| = 2 / pi
        check_area( geode::touching_domain_area( *domain, fset.sampler ),
            100. + 3. * 20. * 2. / M_PI, 0.3, "length in [2, 4], any azimuth" );

        // Small window and long fractures: most of the area is the boundary
        // term, which sampling anchors would estimate badly.
        const auto small_domain =
            geode::build_spatial_domain( square_domain( 2., 50. ) );
        set_sampler( fset, 50., 50., 0., M_PI );
        check_area( geode::touching_domain_area( *small_domain, fset.sampler ),
            4. + 50. * 4. * 2. / M_PI, 0.5,
            "small window, length 50, any azimuth" );

        geode::ObjectSetDefinition< geode::Point2D > pset;
        geode::OpenGeodeStochasticStochasticException::test(
            geode::touching_domain_area( *domain, pset.sampler ) == 100.,
            "[ObservedIntensity] the touching area of points should be the "
            "domain area" );
        geode::Logger::info( "--> SUCCESS!" );
    }

    // lambda counts the fixed points: lambda * |D| = 30 with 10 fixed points,
    // so 20 free points in the domain on average (40 without the correction).
    void test_points_all_objects()
    {
        geode::Logger::info( "TEST - POINTS, DENSITY OF ALL THE OBJECTS" );
        geode::RandomEngine engine;
        engine.set_seed( "@observed-intensity-points@" );

        geode::PoissonProcessBuilder< geode::Point2D > poisson;
        poisson.set_domain( square_domain( 10., 2. ) );
        auto& set_config = poisson.add_set(
            "set_A", 0.3, 30., geode::IntensityScope::all_objects );
        for( const auto i : geode::Range{ 10 } )
        {
            set_config.fixed_objects.emplace_back(
                geode::Point2D{ { 0.5 + i, 5. } } );
        }

        geode::SimulationRunner< geode::Point2D > runner{
            poisson.build_simulation_context()
        };
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 2000;
        sim_config.metropolis_hasting_steps = 100;
        sim_config.burn_in_steps = 1000;
        auto statistic_tracker = runner.run( engine, sim_config );
        geode::TargetStatistics target_stats{ runner.model(),
            poisson.expected_statistics() };
        geode::statistics::validate( statistic_tracker, target_stats );
        geode::Logger::info( "--> SUCCESS!" );
    }

    struct FreeCounts
    {
        double touching{ 0. };
        double not_touching{ 0. };
    };

    FreeCounts mean_free_counts( geode::RandomEngine& engine,
        geode::FractureSimulationRunner& runner,
        const geode::SpatialDomain< 2 >& domain )
    {
        constexpr geode::index_t nb_realizations{ 2000 };
        runner.run( engine, 5000 );
        FreeCounts counts;
        for( const auto realization : geode::Range{ nb_realizations } )
        {
            geode_unused( realization );
            const auto& state = runner.run( engine, 200 );
            const auto& set = state.get_set( state.get_set_uuid( "fset" ) );
            for( const auto id : geode::Range{ set.nb_free_objects() } )
            {
                if( geode::SpatialDomainChecker<
                        geode::Fracture >::touches_domain( domain,
                        set.get_free_object( id ) ) )
                {
                    counts.touching++;
                }
                else
                {
                    counts.not_touching++;
                }
            }
        }
        counts.touching /= nb_realizations;
        counts.not_touching /= nb_realizations;
        return counts;
    }

    // The free fractures touching the domain have the density
    // p20 - F / A, the others keep p20.
    void test_fractures_all_objects()
    {
        geode::Logger::info( "TEST - FRACTURES, DENSITY OF ALL THE OBJECTS" );
        geode::RandomEngine engine;
        engine.set_seed( "@observed-intensity-fractures@" );

        constexpr double p20{ 0.05 };
        const auto domain_cfg = square_domain( 20., 5. );
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( domain_cfg );
        auto& fset = fractures.add_fracture_set(
            "fset", p20, std::nullopt, geode::IntensityScope::all_objects );
        set_sampler( fset, 2., 4., 0., M_PI );
        // 8 observed fractures touching the domain, one crossing its boundary
        for( const auto i : geode::Range{ 7 } )
        {
            fset.fixed_objects.emplace_back(
                geode::Point2D{ { 2. + 2. * i, 5. } },
                geode::Point2D{ { 4. + 2. * i, 8. } } );
        }
        fset.fixed_objects.emplace_back(
            geode::Point2D{ { -2., 10. } }, geode::Point2D{ { 1., 10. } } );

        const auto domain = geode::build_spatial_domain( domain_cfg );
        const auto area = geode::touching_domain_area( *domain, fset.sampler );
        const FreeCounts expected{ p20 * area - 8.,
            p20 * ( domain->extended_n_volume() - area ) };

        geode::FractureSimulationRunner runner{
            fractures.build_simulation_context()
        };
        const auto counts = mean_free_counts( engine, runner, *domain );
        geode::Logger::info( "Free fractures touching the domain: ",
            counts.touching, " (expected: ", expected.touching,
            "), not touching: ", counts.not_touching,
            " (expected: ", expected.not_touching, ")" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( counts.touching - expected.touching )
                < 0.05 * expected.touching,
            "[ObservedIntensity] wrong number of free fractures touching the "
            "domain" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( counts.not_touching - expected.not_touching )
                < 0.05 * expected.not_touching,
            "[ObservedIntensity] wrong number of free fractures not touching "
            "the domain" );
        geode::Logger::info( "--> SUCCESS!" );
    }

    // lambda * |D| = 10 with 10 fixed points: no free point in the domain,
    // free points in the buffer only. A fixed point outside the domain is kept
    // but does not count among the observed ones.
    void test_all_objects_observed()
    {
        geode::Logger::info( "TEST - ALL THE OBJECTS OBSERVED" );
        geode::RandomEngine engine;
        engine.set_seed( "@observed-intensity-all@" );

        geode::PoissonProcessBuilder< geode::Point2D > poisson;
        const auto domain_cfg = square_domain( 10., 2. );
        poisson.set_domain( domain_cfg );
        auto& set_config = poisson.add_set(
            "set_A", 0.1, std::nullopt, geode::IntensityScope::all_objects );
        for( const auto i : geode::Range{ 10 } )
        {
            set_config.fixed_objects.emplace_back(
                geode::Point2D{ { 0.5 + i, 5. } } );
        }
        set_config.fixed_objects.emplace_back( geode::Point2D{ { -1., 5. } } );
        geode::SimulationRunner< geode::Point2D > runner{
            poisson.build_simulation_context()
        };
        const auto domain = geode::build_spatial_domain( domain_cfg );
        geode::index_t nb_in_buffer{ 0 };
        for( const auto realization : geode::Range{ 100 } )
        {
            geode_unused( realization );
            const auto& state = runner.run( engine, 100 );
            const auto& set = state.get_set( state.get_set_uuid( "set_A" ) );
            for( const auto id : geode::Range{ set.nb_free_objects() } )
            {
                const auto& point = set.get_free_object( id );
                geode::OpenGeodeStochasticStochasticException::test(
                    !domain->contains( point ),
                    "[ObservedIntensity] a free point is in the domain whose "
                    "points are all observed: ",
                    point.string() );
                nb_in_buffer++;
            }
        }
        geode::OpenGeodeStochasticStochasticException::test( nb_in_buffer > 0,
            "[ObservedIntensity] no free point in the buffer" );
        geode::Logger::info( "--> SUCCESS!" );
    }

    // lambda * |D| = 5 < 10 fixed points
    void test_underestimated_density()
    {
        geode::Logger::info( "TEST - UNDERESTIMATED DENSITY" );
        geode::PoissonProcessBuilder< geode::Point2D > poisson;
        poisson.set_domain( square_domain( 10., 0. ) );
        auto& set_config = poisson.add_set(
            "set_A", 0.05, std::nullopt, geode::IntensityScope::all_objects );
        for( const auto i : geode::Range{ 10 } )
        {
            set_config.fixed_objects.emplace_back(
                geode::Point2D{ { 0.5 + i, 5. } } );
        }
        try
        {
            const auto context = poisson.build_simulation_context();
        }
        catch( const geode::OpenGeodeException& )
        {
            geode::Logger::info( "--> SUCCESS!" );
            return;
        }
        throw geode::OpenGeodeStochasticStochasticException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[ObservedIntensity] a density lower than the observed objects "
            "should be refused" };
    }
    // NOLINTEND(*-magic-numbers)
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_touching_domain_area();
        test_underestimated_density();
        test_all_objects_observed();
        test_points_all_objects();
        test_fractures_all_objects();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
