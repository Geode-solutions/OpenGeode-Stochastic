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

// Micro-benchmark of the Metropolis-Hastings steps: not a test, nothing is
// checked. For each scenario, the chain is burnt in, then the same number of
// steps is timed several times. Printed: the number of steps per second of the
// fastest repetition (the least disturbed by the other processes of the
// machine, the value to compare), the median one, the acceptance rate and the
// number of objects at the end. Seeds are fixed: two runs do the same work.
//
// Usage: benchmark-metropolis-hastings [scenario name filter]

#include <algorithm>
#include <chrono>
#include <string_view>
#include <vector>

#include <absl/strings/match.h>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/applications/strauss_process.hpp>

namespace
{
    // NOLINTBEGIN(*-magic-numbers)
    constexpr geode::index_t BURN_IN_STEPS{ 200'000 };
    constexpr geode::index_t TIMED_STEPS{ 200'000 };
    constexpr geode::index_t REPETITIONS{ 5 };

    geode::SpatialDomainConfig< 2 > square_domain( double size, double buffer )
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { size, size } };
        domain.buffer_size = buffer;
        return domain;
    }

    void set_fracture_sampler(
        geode::ObjectSetDefinition< geode::Fracture >& fset,
        double min_length,
        double max_length )
    {
        fset.sampler.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.length.min_value = min_length;
        fset.sampler.length.max_value = max_length;
        fset.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.azimuth.min_value = 0.;
        fset.sampler.azimuth.max_value = 180.;
    }

    template < typename ObjectType >
    void run_scenario(
        std::string_view name, geode::SimulationContext< ObjectType > context )
    {
        geode::RandomEngine engine;
        engine.set_seed( name );
        auto& state = *context.object_sets;
        const auto& mh = *context.mh_sampler;
        mh.walk( state, engine, BURN_IN_STEPS );

        std::vector< double > steps_per_second;
        geode::index_t accepted{ 0 };
        for( const auto repetition : geode::Range{ REPETITIONS } )
        {
            geode_unused( repetition );
            const auto start = std::chrono::steady_clock::now();
            for( const auto step : geode::Range{ TIMED_STEPS } )
            {
                geode_unused( step );
                if( mh.step( state, engine ).decision
                    == geode::MH_DECISION::accepted )
                {
                    accepted++;
                }
            }
            const std::chrono::duration< double > elapsed =
                std::chrono::steady_clock::now() - start;
            steps_per_second.push_back( TIMED_STEPS / elapsed.count() );
        }
        std::sort( steps_per_second.begin(), steps_per_second.end() );
        const auto acceptance =
            static_cast< double >( accepted ) / ( REPETITIONS * TIMED_STEPS );
        const auto& final_state = state;
        geode::index_t nb_free{ 0 };
        for( const auto& set_name : final_state.set_names() )
        {
            nb_free +=
                final_state.get_set( final_state.get_set_uuid( set_name ) )
                    .nb_free_objects();
        }
        geode::Logger::info( name, ": ",
            static_cast< int >( steps_per_second.back() ), " steps/s (median ",
            static_cast< int >( steps_per_second[REPETITIONS / 2] ),
            "), acceptance ", static_cast< int >( 100. * acceptance ), " %, ",
            nb_free, " free objects, ", state.nb_objects() - nb_free,
            " fixed" );
    }

    // Strauss 2D: about 500 points in the domain
    geode::SimulationContext< geode::Point2D > strauss_context()
    {
        geode::StraussProcessBuilder< geode::Point2D > strauss;
        strauss.set_domain( square_domain( 40., 2. ) );
        auto& set = strauss.add_set( "points", 0.5 );
        geode_unused( set );
        strauss.add_interaction( { "points" }, 0.5, 1., std::nullopt );
        return strauss.build_simulation_context();
    }

    // Two fracture sets, with an interaction between intersecting fractures
    // of different sets.
    geode::FractureSimulationContext x_node_context()
    {
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( square_domain( 20., 3. ) );
        for( const auto* set_name : { "fset_01", "fset_02" } )
        {
            auto& fset = fractures.add_fracture_set( set_name, 0.3 );
            set_fracture_sampler( fset, 1., 3. );
        }
        fractures.add_x_node_interaction( { "fset_01", "fset_02" }, 0.5 );
        return fractures.build_simulation_context();
    }

    // Small fractures that cannot intersect, with large fixed fractures
    // parallel to the diagonal of the domain, crossing the whole extended
    // domain: the bounding box of each one covers a large part of it.
    geode::FractureSimulationContext large_fixed_context(
        geode::index_t nb_large )
    {
        constexpr double size{ 20. };
        constexpr double buffer{ 1. };
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( square_domain( size, buffer ) );
        auto& fset = fractures.add_fracture_set( "fset", 0.3 );
        set_fracture_sampler( fset, 0.5, 1. );
        // line x - y = offset, clipped to the extended domain
        constexpr double low{ -buffer };
        constexpr double high{ size + buffer };
        for( const auto line : geode::Range{ nb_large } )
        {
            const auto offset =
                ( high - low ) * ( 2. * ( line + 0.5 ) / nb_large - 1. );
            if( offset >= 0. )
            {
                fset.fixed_objects.emplace_back(
                    geode::Point2D{ { low + offset, low } },
                    geode::Point2D{ { high, high - offset } } );
            }
            else
            {
                fset.fixed_objects.emplace_back(
                    geode::Point2D{ { low, low - offset } },
                    geode::Point2D{ { high + offset, high } } );
            }
        }
        fractures.add_minimal_spacing( "fset", 0. );
        return fractures.build_simulation_context();
    }
    // NOLINTEND(*-magic-numbers)
} // namespace

int main( int argc, char* argv[] )
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        const std::string_view filter{ argc > 1 ? argv[1] : "" };
        const auto selected = [filter]( std::string_view name ) {
            return absl::StrContains( name, filter );
        };
        if( selected( "strauss" ) )
        {
            run_scenario( "strauss", strauss_context() );
        }
        if( selected( "x_node" ) )
        {
            run_scenario( "x_node", x_node_context() );
        }
        for( const geode::index_t nb_large : { 0, 10, 30, 100, 300 } )
        {
            const auto name = absl::StrCat( "large_fixed_", nb_large );
            if( selected( name ) )
            {
                run_scenario( name, large_fixed_context( nb_large ) );
            }
        }
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
