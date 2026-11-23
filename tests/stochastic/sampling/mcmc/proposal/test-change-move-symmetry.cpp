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
#include <memory>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/applications/poisson_process.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/moves.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

// Symmetry test of the point change move: a Poisson process simulated with
// change moves only (fixed number of points, no interaction) has a uniform
// stationary law. A change kernel that is not symmetric near the boundary of
// the extended domain makes the density there differ from the density in the
// center.
namespace
{
    constexpr double DOMAIN_SIZE{ 20. };
    // move step = MOVE_RATIO * DOMAIN_SIZE, equal to the strip width: the
    // whole strip is within one step of the boundary
    constexpr double MOVE_RATIO{ 0.1 };
    constexpr double STRIP_WIDTH{ MOVE_RATIO * DOMAIN_SIZE };
    constexpr double CENTER_HALF_SIZE{ 5. };
    constexpr geode::index_t NB_POINTS{ 200 };
    // relative difference between strip and center densities accepted for a
    // uniform law (statistical fluctuations); a truncated kernel gives about
    // -20%
    constexpr double TOLERANCE{ 0.04 };

    geode::SimulationContext< geode::Point2D > build_change_only_context(
        geode::RandomEngine& engine )
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        domain.buffer_size = 0.;

        geode::PoissonProcessBuilder< geode::Point2D > poisson;
        poisson.set_domain( domain );
        auto& set = poisson.add_set( "points", 1. );
        set.sampler.move_ratio = MOVE_RATIO;
        auto context = poisson.build_simulation_context();

        const auto set_id = context.object_sets->get_set_uuid( "points" );
        const auto& sampler = *context.set_samplers.front();
        for( const auto point : geode::Range{ NB_POINTS } )
        {
            geode_unused( point );
            context.object_sets->add_object(
                sampler.sample( engine ), set_id, false );
        }
        auto kernel =
            std::make_unique< geode::ProposalKernel< geode::Point2D > >();
        kernel->add_move(
            set_id, std::make_unique< geode::ChangeMove< geode::Point2D > >(
                        sampler, 1. ) );
        context.mh_sampler =
            std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                *context.model, std::move( kernel ) );
        return context;
    }

    void test_point_change_move_symmetry()
    {
        geode::Logger::info( "TEST - POINT CHANGE MOVE SYMMETRY" );
        geode::RandomEngine engine;
        engine.set_seed( "@point-change-move-symmetry@" );

        // NOLINTBEGIN(*-magic-numbers)
        constexpr geode::index_t burn_in_steps{ 20000 };
        constexpr geode::index_t realizations{ 1000 };
        constexpr geode::index_t steps{ 1000 };
        // NOLINTEND(*-magic-numbers)

        geode::SimulationRunner< geode::Point2D > runner{
            build_change_only_context( engine )
        };
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
            const auto objects =
                state.get_objects_in_set( state.get_set_uuid( "points" ) );
            geode::OpenGeodeStochasticStochasticException::test(
                objects.size() == NB_POINTS,
                "[ChangeMoveSymmetry] The number of points should not change" );
            for( const auto& object_id : objects )
            {
                const auto& point = state.get_object( object_id );
                const auto x = point.value( 0 );
                const auto y = point.value( 1 );
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
        const auto inner_size = DOMAIN_SIZE - 2. * STRIP_WIDTH;
        const auto strip_area =
            DOMAIN_SIZE * DOMAIN_SIZE - inner_size * inner_size;
        const auto center_area = 4. * CENTER_HALF_SIZE * CENTER_HALF_SIZE;
        const auto strip_density = nb_strip / strip_area;
        const auto center_density = nb_center / center_area;
        const auto relative_difference =
            ( strip_density - center_density ) / center_density;
        geode::Logger::info( "strip density=", strip_density / realizations,
            " center density=", center_density / realizations,
            " relative difference=", relative_difference );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( relative_difference ) < TOLERANCE,
            "[ChangeMoveSymmetry] Non uniform point density near the domain "
            "boundary (relative difference = ",
            relative_difference, ")" );
        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        test_point_change_move_symmetry();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
