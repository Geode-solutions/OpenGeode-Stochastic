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

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/applications/strauss_process.hpp>

// Energy consistency test: along a random sequence of births, deaths and
// changes, the variation of the total energy must equal the energy deltas
// used by the Metropolis-Hastings sampler, for objects anywhere in the
// extended domain (domain + buffer).
namespace
{
    constexpr geode::index_t NB_OPERATIONS{ 1000 };

    geode::index_t random_index(
        geode::RandomEngine& engine, geode::index_t size )
    {
        geode::UniformClosedOpen< double > law;
        law.min_value = 0.;
        law.max_value = static_cast< double >( size );
        return static_cast< geode::index_t >(
            std::floor( engine.sample_uniform( law ) ) );
    }

    void check_delta(
        double before, double after, double delta, std::string_view operation )
    {
        const auto error = std::fabs( after - before - delta );
        geode::OpenGeodeStochasticStochasticException::test(
            error < 1e-9 * ( 1. + std::fabs( after ) ), "[EnergyConsistency] ",
            operation, " delta = ", delta,
            " but total energy variation = ", after - before );
    }

    template < typename ObjectType >
    void check_energy_consistency( geode::RandomEngine& engine,
        geode::SimulationContext< ObjectType >& context,
        const std::vector< std::string >& set_names )
    {
        auto& state = *context.object_sets;
        const auto& energy = context.model->energy();
        for( const auto operation : geode::Range{ NB_OPERATIONS } )
        {
            geode_unused( operation );
            const auto set_index = random_index( engine, set_names.size() );
            const auto set_id = state.get_set_uuid( set_names[set_index] );
            const auto& sampler = *context.set_samplers[set_index];
            const auto before = energy.total_log_energy( state );
            const auto objects = state.get_objects_in_set( set_id );

            // NOLINTNEXTLINE(*-magic-numbers)
            if( objects.empty() || engine.sample_bernoulli( 0.5 ) )
            {
                auto object = sampler.sample( engine );
                const auto delta = energy.delta_log_add(
                    state, geode::ObjectRef< ObjectType >{ object, set_id } );
                state.add_object( std::move( object ), set_id, false );
                check_delta(
                    before, energy.total_log_energy( state ), delta, "add" );
                continue;
            }
            const auto& object_id =
                objects[random_index( engine, objects.size() )];
            if( engine.sample_bernoulli( 0.5 ) )
            {
                const auto delta = energy.delta_log_remove( state, object_id );
                state.remove_free_object( object_id );
                check_delta(
                    before, energy.total_log_energy( state ), delta, "remove" );
                continue;
            }
            auto object =
                sampler.change( state.get_object( object_id ), engine );
            const auto delta = energy.delta_log_change( state, object_id,
                geode::ObjectRef< ObjectType >{ object, set_id } );
            state.update_free_object( object_id, std::move( object ) );
            check_delta(
                before, energy.total_log_energy( state ), delta, "change" );
        }
    }

    geode::SpatialDomainConfig< 2 > domain_config()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { 5., 5. } };
        domain.buffer_size = 2.;
        // NOLINTEND(*-magic-numbers)
        return domain;
    }

    void test_point_process()
    {
        geode::Logger::info( "TEST - ENERGY CONSISTENCY POINT PROCESS" );
        geode::RandomEngine engine;
        engine.set_seed( "@energy-consistency-points@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::StraussProcessBuilder< geode::Point2D > strauss;
        strauss.set_domain( domain_config() );
        auto& set_a = strauss.add_set( "A", 2. );
        geode_unused( set_a );
        auto& set_b = strauss.add_set( "B", 1. );
        geode_unused( set_b );
        strauss.add_interaction( { "A", "B" }, 0.5, 1., std::nullopt, true );
        strauss.add_interaction( { "A", "B" }, 0.8, 0.5, std::nullopt, false );
        // NOLINTEND(*-magic-numbers)

        auto context = strauss.build_simulation_context();
        check_energy_consistency( engine, context, { "A", "B" } );
        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_segment_process()
    {
        geode::Logger::info( "TEST - ENERGY CONSISTENCY SEGMENT PROCESS" );
        geode::RandomEngine engine;
        engine.set_seed( "@energy-consistency-segments@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( domain_config() );
        for( const auto& name : { "A", "B" } )
        {
            auto& fset = fractures.add_fracture_set( name, 1. );
            fset.sampler.length.distribution_type =
                geode::UniformClosed< double >::distribution_type_static();
            fset.sampler.length.min_value = 0.5;
            fset.sampler.length.max_value = 2.;
            fset.sampler.azimuth.distribution_type =
                geode::UniformClosed< double >::distribution_type_static();
            fset.sampler.azimuth.min_value = 0.;
            fset.sampler.azimuth.max_value = 180.;
            fractures.add_intensity( name, 1.5 );
        }
        fractures.add_x_node_interaction( { "A", "B" }, 0.5 );
        // NOLINTEND(*-magic-numbers)

        auto context = fractures.build_simulation_context();
        check_energy_consistency( engine, context, { "A", "B" } );
        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_point_process();
        test_segment_process();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
