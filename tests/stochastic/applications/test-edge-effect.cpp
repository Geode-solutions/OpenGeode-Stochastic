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

#include <algorithm>

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/applications/strauss_process.hpp>

#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

// Edge effect test: for a stationary process, the density of objects anchored
// near the domain boundary must equal the density in the domain center.
// Missing (no buffer) or wrongly simulated (buffer objects not interacting)
// neighbors outside the domain make the process inhomogeneous in the domain.
namespace
{
    constexpr double DOMAIN_SIZE{ 20. };
    // width of the boundary strip of the domain
    constexpr double STRIP_WIDTH{ 2. };
    // half size of the central window
    constexpr double CENTER_HALF_SIZE{ 5. };

    geode::SpatialDomainConfig< 2 > domain_config( double buffer )
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        domain.buffer_size = buffer;
        return domain;
    }

    const geode::Point2D& anchor( const geode::Point2D& point )
    {
        return point;
    }

    const geode::Point2D& anchor( const geode::OwnerSegment2D& segment )
    {
        return segment.vertices()[0];
    }

    struct EdgeStatistics
    {
        double strip_density;
        double center_density;

        [[nodiscard]] double relative_difference() const
        {
            return ( strip_density - center_density ) / center_density;
        }
    };

    template < typename ObjectType >
    EdgeStatistics measure_edge_statistics( geode::RandomEngine& engine,
        geode::SimulationContext< ObjectType >&& context,
        const std::vector< std::string >& set_names )
    {
        // NOLINTBEGIN(*-magic-numbers)
        constexpr geode::index_t burn_in_steps{ 20000 };
        constexpr geode::index_t realizations{ 1000 };
        constexpr geode::index_t steps{ 1000 };
        // NOLINTEND(*-magic-numbers)

        geode::SimulationRunner< ObjectType > runner{ std::move( context ) };
        const auto& burn_state = runner.run( engine, burn_in_steps );
        geode_unused( burn_state );

        constexpr double center_min{ DOMAIN_SIZE / 2. - CENTER_HALF_SIZE };
        constexpr double center_max{ DOMAIN_SIZE / 2. + CENTER_HALF_SIZE };
        double nb_strip{ 0. };
        double nb_center{ 0. };
        for( const auto realization : geode::Range{ realizations } )
        {
            geode_unused( realization );
            const auto& state = runner.run( engine, steps );
            for( const auto& set_name : set_names )
            {
                for( const auto& object_id :
                    state.get_objects_in_set( state.get_set_uuid( set_name ) ) )
                {
                    const auto& point = anchor( state.get_object( object_id ) );
                    const auto x = point.value( 0 );
                    const auto y = point.value( 1 );
                    if( x < 0. || y < 0. || x > DOMAIN_SIZE || y > DOMAIN_SIZE )
                    {
                        continue;
                    }
                    const auto boundary_distance = std::min< double >(
                        { x, y, DOMAIN_SIZE - x, DOMAIN_SIZE - y } );
                    if( boundary_distance < STRIP_WIDTH )
                    {
                        nb_strip += 1.;
                    }
                    if( x > center_min && x < center_max && y > center_min
                        && y < center_max )
                    {
                        nb_center += 1.;
                    }
                }
            }
        }
        const auto inner_size = DOMAIN_SIZE - 2. * STRIP_WIDTH;
        const auto strip_area =
            DOMAIN_SIZE * DOMAIN_SIZE - inner_size * inner_size;
        const auto center_area = 4. * CENTER_HALF_SIZE * CENTER_HALF_SIZE;
        return { nb_strip / ( realizations * strip_area ),
            nb_center / ( realizations * center_area ) };
    }

    // relative difference between boundary and center densities accepted for
    // an unbiased process (statistical fluctuations)
    constexpr double TOLERANCE{ 0.04 };

    void check_statistics( double buffer,
        const EdgeStatistics& stats,
        bool expect_edge_effect,
        std::string_view process )
    {
        geode::Logger::info( process, " - buffer=", buffer,
            " boundary density=", stats.strip_density,
            " center density=", stats.center_density,
            " relative difference=", stats.relative_difference() );
        const auto has_edge_effect =
            std::fabs( stats.relative_difference() ) > TOLERANCE;
        geode::OpenGeodeStochasticStochasticException::test(
            has_edge_effect == expect_edge_effect, "[EdgeEffect] ", process,
            " with buffer = ", buffer,
            expect_edge_effect ? ": edge effect expected but not detected"
                               : ": biased process near the domain boundary",
            " (relative difference = ", stats.relative_difference(), ")" );
    }

    void set_uniform_sampler(
        geode::ObjectSetDefinition< geode::Fracture >& fset )
    {
        // NOLINTBEGIN(*-magic-numbers)
        fset.sampler.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.length.min_value = 0.5;
        fset.sampler.length.max_value = 1.;
        fset.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.azimuth.min_value = 0.;
        fset.sampler.azimuth.max_value = 180.;
        // NOLINTEND(*-magic-numbers)
    }

    void test_hard_core_points()
    {
        geode::Logger::info( "TEST - EDGE EFFECT HARD-CORE POINT PROCESS" );
        geode::RandomEngine engine;
        engine.set_seed( "@edge-effect-hard-core-points@" );

        // NOLINTBEGIN(*-magic-numbers)
        // without buffer, points near the boundary miss their outside
        // neighbors: the edge effect must be detected (test sensitivity)
        for( const double buffer : { 0., 2., 4. } )
        {
            geode::StraussProcessBuilder< geode::Point2D > strauss;
            strauss.set_domain( domain_config( buffer ) );
            auto& set = strauss.add_set( "points", 1. );
            geode_unused( set );
            strauss.add_interaction( { "points" }, 0., 1., std::nullopt );
            check_statistics( buffer,
                measure_edge_statistics(
                    engine, strauss.build_simulation_context(), { "points" } ),
                buffer == 0., "hard-core points" );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_hard_core_segments()
    {
        geode::Logger::info( "TEST - EDGE EFFECT HARD-CORE SEGMENT PROCESS" );
        geode::RandomEngine engine;
        engine.set_seed( "@edge-effect-hard-core-segments@" );

        // NOLINTBEGIN(*-magic-numbers)
        for( const double buffer : { 3., 5. } )
        {
            geode::FractureProcessBuilder fractures;
            fractures.set_domain( domain_config( buffer ) );
            auto& fset = fractures.add_fracture_set( "fset", 1. );
            set_uniform_sampler( fset );
            fractures.add_minimal_spacing( "fset", 0.5 );
            check_statistics( buffer,
                measure_edge_statistics(
                    engine, fractures.build_simulation_context(), { "fset" } ),
                false, "hard-core segments" );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_inhibited_x_nodes()
    {
        geode::Logger::info( "TEST - EDGE EFFECT INHIBITED X-NODES" );
        geode::RandomEngine engine;
        engine.set_seed( "@edge-effect-x-nodes@" );

        // NOLINTBEGIN(*-magic-numbers)
        for( const double buffer : { 3., 5. } )
        {
            geode::FractureProcessBuilder fractures;
            fractures.set_domain( domain_config( buffer ) );
            auto& fset_a = fractures.add_fracture_set( "fset_A", 1. );
            auto& fset_b = fractures.add_fracture_set( "fset_B", 1. );
            set_uniform_sampler( fset_a );
            set_uniform_sampler( fset_b );
            fractures.add_x_node_interaction( { "fset_A", "fset_B" }, 0. );
            check_statistics( buffer,
                measure_edge_statistics( engine,
                    fractures.build_simulation_context(),
                    { "fset_A", "fset_B" } ),
                false, "inhibited X-nodes" );
        }
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_hard_core_points();
        test_hard_core_segments();
        test_inhibited_x_nodes();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
