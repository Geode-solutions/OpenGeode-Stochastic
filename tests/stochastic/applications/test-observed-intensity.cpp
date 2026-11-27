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

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/distance.hpp>

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
        const auto& burn_in_state = runner.run( engine, 5000 );
        geode_unused( burn_in_state );
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

    struct LargeFracturesRun
    {
        double mean_free{ 0. };
        geode::index_t nb_crossings{ 0 };
    };

    // Large observed fractures, longer than the buffer and crossing each
    // other, fixed in the set of the small simulated ones.
    LargeFracturesRun run_with_large_fractures(
        std::string_view seed, bool forbid_intersections )
    {
        geode::RandomEngine engine;
        engine.set_seed( seed );
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( square_domain( 20., 5. ) );
        auto& fset = fractures.add_fracture_set( "fset", 0.05 );
        // the length law of the simulated fractures is cut at L_obs = 4
        set_sampler( fset, 1., 4., 0., M_PI );
        std::vector< geode::Fracture > large;
        large.emplace_back(
            geode::Point2D{ { -10., 7. } }, geode::Point2D{ { 30., 7. } } );
        large.emplace_back(
            geode::Point2D{ { -10., 14. } }, geode::Point2D{ { 30., 14. } } );
        large.emplace_back(
            geode::Point2D{ { -8., -8. } }, geode::Point2D{ { 28., 28. } } );
        fset.fixed_objects = large;
        if( forbid_intersections )
        {
            fractures.add_minimal_spacing( "fset", 0. );
        }
        geode::FractureSimulationRunner runner{
            fractures.build_simulation_context()
        };

        constexpr geode::index_t nb_realizations{ 50 };
        const auto& burn_in_state = runner.run( engine, 5000 );
        geode_unused( burn_in_state );
        LargeFracturesRun result;
        for( const auto realization : geode::Range{ nb_realizations } )
        {
            geode_unused( realization );
            const auto& state = runner.run( engine, 200 );
            const auto& set = state.get_set( state.get_set_uuid( "fset" ) );
            geode::OpenGeodeStochasticStochasticException::test(
                set.nb_fixed_objects() == large.size(),
                "[ObservedIntensity] the large fractures changed" );
            for( const auto fixed_id : geode::Range{ large.size() } )
            {
                const auto& fixed = set.get_fixed_object( fixed_id );
                geode::OpenGeodeStochasticStochasticException::test(
                    fixed.vertices() == large[fixed_id].vertices(),
                    "[ObservedIntensity] a large fracture moved" );
            }
            result.mean_free += set.nb_free_objects();
            for( const auto id : geode::Range{ set.nb_free_objects() } )
            {
                const geode::Segment2D small{ set.get_free_object( id ) };
                for( const auto& fixed : large )
                {
                    if( std::get< 0 >( geode::segment_segment_distance(
                            small, geode::Segment2D{ fixed } ) )
                        < geode::GLOBAL_EPSILON )
                    {
                        result.nb_crossings++;
                    }
                }
            }
        }
        result.mean_free /= nb_realizations;
        return result;
    }

    // The large fractures are observed and fixed in the set of the small
    // simulated ones, which see them through the interaction terms of the set
    // (here no intersection), even if they are longer than the buffer.
    void test_large_fractures_interaction()
    {
        geode::Logger::info( "TEST - LARGE OBSERVED FRACTURES, INTERACTION" );
        const auto free_run =
            run_with_large_fractures( "@large-fractures-free@", false );
        geode::Logger::info( "Without interaction: ", free_run.mean_free,
            " simulated fractures, ", free_run.nb_crossings,
            " crossing a large one" );
        geode::OpenGeodeStochasticStochasticException::test(
            free_run.nb_crossings > 0, "[ObservedIntensity] no small fracture "
                                       "crosses a large one without "
                                       "interaction" );

        const auto spacing_run =
            run_with_large_fractures( "@large-fractures-spacing@", true );
        geode::Logger::info( "Without intersection: ", spacing_run.mean_free,
            " simulated fractures, ", spacing_run.nb_crossings,
            " crossing a large one" );
        geode::OpenGeodeStochasticStochasticException::test(
            spacing_run.nb_crossings == 0,
            "[ObservedIntensity] small fractures cross the large ones despite "
            "the minimal spacing" );
        geode::OpenGeodeStochasticStochasticException::test(
            spacing_run.mean_free > 0.,
            "[ObservedIntensity] no small fracture simulated" );
        geode::Logger::info( "--> SUCCESS!" );
    }

    void check_free_counts( const FreeCounts& counts,
        const FreeCounts& expected,
        std::string_view label )
    {
        geode::Logger::info( label,
            ": free fractures touching the domain: ", counts.touching,
            " (expected: ", expected.touching,
            "), not touching: ", counts.not_touching,
            " (expected: ", expected.not_touching, ")" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( counts.touching - expected.touching )
                < 0.05 * expected.touching,
            "[ObservedIntensity] ", label,
            ": wrong number of free fractures touching the domain" );
        geode::OpenGeodeStochasticStochasticException::test(
            std::abs( counts.not_touching - expected.not_touching )
                < 0.05 * expected.not_touching,
            "[ObservedIntensity] ", label,
            ": wrong number of free fractures not touching the domain" );
    }

    // One fracture network, p20 = 0.1 with lengths uniform in [1, 8],
    // observed in three ways. Large fractures are longer than L_obs = 4.
    // 1. A sample of all sizes: length law of all the fractures.
    // 2. All the large fractures, no small one: length law cut at L_obs,
    // p20 of the small fractures (0.1 * 3 / 7); no large fracture counts
    // as observed (F = 0), as with free_objects.
    // 3. All the large fractures and a sample of the small ones: only the
    // small observed fractures count in F.
    void test_three_observation_cases()
    {
        geode::Logger::info( "TEST - THREE OBSERVATION CASES" );
        constexpr double p20{ 0.1 };
        constexpr double l_obs{ 4. };
        const auto small_p20 = p20 * ( l_obs - 1. ) / ( 8. - 1. );
        const auto domain_cfg = square_domain( 20., 8. );
        const auto domain = geode::build_spatial_domain( domain_cfg );
        const auto extended_area = domain->extended_n_volume();

        // observed small fractures (length 3) and large ones (lengths 6,
        // and 30, longer than the buffer), all touching the domain
        std::vector< geode::Fracture > small_observed;
        for( const auto i : geode::Range{ 4 } )
        {
            small_observed.emplace_back( geode::Point2D{ { 2. + 4. * i, 3. } },
                geode::Point2D{ { 2. + 4. * i, 6. } } );
        }
        std::vector< geode::Fracture > large_observed;
        large_observed.emplace_back(
            geode::Point2D{ { 3., 10. } }, geode::Point2D{ { 9., 10. } } );
        large_observed.emplace_back(
            geode::Point2D{ { -5., 15. } }, geode::Point2D{ { 25., 15. } } );

        const auto run_case = [&]( std::string_view seed, double density,
                                  double max_length,
                                  const std::vector< geode::Fracture >& fixed,
                                  std::string_view label, double nb_observed ) {
            geode::RandomEngine engine;
            engine.set_seed( seed );
            geode::FractureProcessBuilder fractures;
            fractures.set_domain( domain_cfg );
            auto& fset = fractures.add_fracture_set( "fset", density,
                std::nullopt, geode::IntensityScope::all_objects );
            set_sampler( fset, 1., max_length, 0., M_PI );
            fset.fixed_objects = fixed;
            const auto area =
                geode::touching_domain_area( *domain, fset.sampler );
            geode::FractureSimulationRunner runner{
                fractures.build_simulation_context()
            };
            check_free_counts( mean_free_counts( engine, runner, *domain ),
                { density * area - nb_observed,
                    density * ( extended_area - area ) },
                label );
        };

        // 1. the sample mixes small and large fractures
        auto sample = small_observed;
        sample.push_back( large_observed.front() );
        run_case( "@observation-case-1@", p20, 8., sample, "1. sample", 5. );

        // 2. all the large fractures, no small one
        run_case( "@observation-case-2@", small_p20, l_obs, large_observed,
            "2. all large", 0. );

        // 3. all the large fractures and a sample of the small ones
        auto both = large_observed;
        both.insert( both.end(), small_observed.begin(), small_observed.end() );
        run_case( "@observation-case-3@", small_p20, l_obs, both,
            "3. all large and small sample", 4. );
        geode::Logger::info( "--> SUCCESS!" );
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
        test_three_observation_cases();
        test_large_fractures_interaction();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
