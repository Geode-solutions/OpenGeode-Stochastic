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

#include <geode/stochastic/inference/statistics_tools.hpp>

namespace
{
    constexpr double DOMAIN_SIZE{ 20. };
    constexpr double DOMAIN_AREA{ DOMAIN_SIZE * DOMAIN_SIZE };
    constexpr double MIN_LENGTH{ 0.5 };
    constexpr double MAX_LENGTH{ 1. };

    geode::SpatialDomainConfig< 2 > fracture_domain()
    {
        geode::SpatialDomainConfig< 2 > domain;
        domain.min_point = geode::Point2D{ { 0., 0. } };
        domain.max_point = geode::Point2D{ { DOMAIN_SIZE, DOMAIN_SIZE } };
        // buffer larger than the longest fracture: no edge effect
        domain.buffer_size = 2. * MAX_LENGTH;
        return domain;
    }

    void set_uniform_sampler(
        geode::ObjectSetDefinition< geode::Fracture >& fset )
    {
        // NOLINTBEGIN(*-magic-numbers)
        fset.sampler.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.length.min_value = MIN_LENGTH;
        fset.sampler.length.max_value = MAX_LENGTH;

        fset.sampler.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        fset.sampler.azimuth.min_value = 0.;
        fset.sampler.azimuth.max_value = 180.;
        // NOLINTEND(*-magic-numbers)
    }

    // For L ~ U[MIN_LENGTH, MAX_LENGTH] and an unnormalized density
    // p20 * p21^L (characteristic length = 1), returns the expected number
    // and total length of fractures in the domain.
    std::pair< double, double > expected_statistics( double p20, double p21 )
    {
        const auto range = MAX_LENGTH - MIN_LENGTH;
        if( std::fabs( p21 - 1. ) < geode::GLOBAL_EPSILON )
        {
            const auto nb = p20 * DOMAIN_AREA;
            return { nb, nb * ( MIN_LENGTH + MAX_LENGTH ) / 2. };
        }
        const auto log_p21 = std::log( p21 );
        const auto primitive_weight = [log_p21]( double length ) {
            return std::exp( log_p21 * length ) / log_p21;
        };
        const auto primitive_length = [log_p21]( double length ) {
            return std::exp( log_p21 * length )
                   * ( length / log_p21 - 1. / ( log_p21 * log_p21 ) );
        };
        const auto mean_weight =
            ( primitive_weight( MAX_LENGTH ) - primitive_weight( MIN_LENGTH ) )
            / range;
        const auto mean_weighted_length =
            ( primitive_length( MAX_LENGTH ) - primitive_length( MIN_LENGTH ) )
            / range;
        return { p20 * DOMAIN_AREA * mean_weight,
            p20 * DOMAIN_AREA * mean_weighted_length };
    }

    void run_and_validate( geode::RandomEngine& engine,
        const geode::FractureProcessBuilder& fractures,
        std::string_view output_name )
    {
        geode::FractureSimulationRunner runner{
            fractures.build_simulation_context()
        };

        // NOLINTBEGIN(*-magic-numbers)
        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 1000;
        sim_config.metropolis_hasting_steps = 100;
        sim_config.burn_in_steps = 1000;
        // NOLINTEND(*-magic-numbers)

        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder =
            absl::StrCat( printer_config.output_folder, "/", output_name );
        sim_config.printer = printer_config;

        auto statistic_tracker = runner.run( engine, sim_config );
        geode::TargetStatistics target_stats{ runner.model(),
            fractures.expected_statistics() };
        geode::statistics::validate( statistic_tracker, target_stats );
    }

    void test_fracture_density()
    {
        geode::Logger::info( "TEST - FRACTURE SET P20" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-fracture-p20@" );

        // NOLINTBEGIN(*-magic-numbers)
        constexpr double p20{ 0.05 };
        const auto [nb, length] = expected_statistics( p20, 1. );

        geode::FractureProcessBuilder fractures;
        fractures.set_domain( fracture_domain() );
        auto& fset = fractures.add_fracture_set( "fset_A", p20, nb );
        set_uniform_sampler( fset );
        // p21 = 1: no effect on the model, used to monitor the total length
        fractures.add_intensity( "fset_A", 1., length );
        // NOLINTEND(*-magic-numbers)

        run_and_validate( engine, fractures, "sim_fracture_p20_test" );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_fracture_intensity()
    {
        geode::Logger::info( "TEST - FRACTURE SET P20 + P21" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-fracture-p21@" );

        // NOLINTBEGIN(*-magic-numbers)
        constexpr double p20{ 0.02 };
        constexpr double p21{ 4. };
        const auto [nb, length] = expected_statistics( p20, p21 );

        // change moves: translation / rotation / stretch weights. Each move
        // alone must preserve the target (stretch changes the lengths).
        const std::array< std::array< double, 3 >, 4 > move_ratios{
            { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. }, { 1., 1., 1. } }
        };
        for( const auto& ratios : move_ratios )
        {
            geode::FractureProcessBuilder fractures;
            fractures.set_domain( fracture_domain() );
            auto& fset = fractures.add_fracture_set( "fset_A", p20, nb );
            set_uniform_sampler( fset );
            fset.sampler.translation_ratio = ratios[0];
            fset.sampler.rotation_ratio = ratios[1];
            fset.sampler.stretch_ratio = ratios[2];
            // change moves dominate birth / death moves
            fset.dynamics.change_ratio = 10.;
            fractures.add_intensity( "fset_A", p21, length );

            geode::Logger::info( "change moves (translation / rotation / "
                                 "stretch): ",
                ratios[0], " / ", ratios[1], " / ", ratios[2] );
            run_and_validate( engine, fractures, "sim_fracture_p21_test" );
        }
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_fracture_minimal_spacing()
    {
        geode::Logger::info( "TEST - FRACTURE SET MINIMAL SPACING" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-fracture-spacing@" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::FractureProcessBuilder fractures;
        fractures.set_domain( fracture_domain() );
        auto& fset = fractures.add_fracture_set( "fset_A", 0.1 );
        set_uniform_sampler( fset );
        // hard-core: no pair of fractures closer than the spacing
        fractures.add_minimal_spacing( "fset_A", 1., 0. );
        // NOLINTEND(*-magic-numbers)

        run_and_validate( engine, fractures, "sim_fracture_spacing_test" );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_observed_fractures()
    {
        geode::Logger::info( "TEST - FRACTURE SET WITH OBSERVED FRACTURES" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-fracture-observed@" );

        // NOLINTBEGIN(*-magic-numbers)
        const std::array< geode::Fracture, 2 > observed{
            geode::Fracture{
                geode::Point2D{ { 0., 15. } }, geode::Point2D{ { 15., 15. } } },
            geode::Fracture{
                geode::Point2D{ { 1., 11. } }, geode::Point2D{ { 11., 20. } } }
        };

        geode::FractureProcessBuilder fractures;
        fractures.set_domain( fracture_domain() );
        auto& fset = fractures.add_fracture_set( "fset_A", 0.05 );
        set_uniform_sampler( fset );
        for( const auto& fracture : observed )
        {
            fset.fixed_objects.push_back( fracture );
        }
        // NOLINTEND(*-magic-numbers)

        geode::FractureSimulationRunner runner{
            fractures.build_simulation_context()
        };
        const auto& state = runner.run( engine, 1000 );

        const auto& fset_state =
            state.get_set( state.get_set_uuid( "fset_A" ) );
        geode::OpenGeodeStochasticStochasticException::test(
            fset_state.nb_fixed_objects() == observed.size(),
            "[Fractures] wrong number of observed fractures: ",
            fset_state.nb_fixed_objects() );
        for( const auto id : geode::Range{ observed.size() } )
        {
            const auto& vertices = fset_state.get_fixed_object( id ).vertices();
            const auto& expected = observed[id].vertices();
            geode::OpenGeodeStochasticStochasticException::test(
                vertices[0] == expected[0] && vertices[1] == expected[1],
                "[Fractures] observed fracture ", id,
                " has been modified by the simulation" );
        }

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_two_fracture_sets()
    {
        geode::Logger::info( "TEST - TWO FRACTURE SETS (X-node interaction)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-two-fracture-sets@" );

        // NOLINTBEGIN(*-magic-numbers)
        constexpr double p20_01{ 0.05 };
        constexpr double p20_02{ 0.03 };
        for( const double beta : { 1., 0. } )
        {
            geode::FractureProcessBuilder fractures;
            fractures.set_domain( fracture_domain() );

            // beta = 1: no interaction, each set is an independent process
            const auto independent = beta == 1.;
            const auto expected_count =
                [independent]( double p20 ) -> std::optional< double > {
                if( !independent )
                {
                    return std::nullopt;
                }
                return expected_statistics( p20, 1. ).first;
            };
            auto& fset_01 = fractures.add_fracture_set(
                "fset_01", p20_01, expected_count( p20_01 ) );
            set_uniform_sampler( fset_01 );

            auto& fset_02 = fractures.add_fracture_set(
                "fset_02", p20_02, expected_count( p20_02 ) );
            set_uniform_sampler( fset_02 );
            fset_02.sampler.azimuth.distribution_type =
                geode::VonMises::distribution_type_static();
            fset_02.sampler.azimuth.mean = 60.;
            fset_02.sampler.azimuth.kappa = 1.;

            // beta = 0: X-nodes are inhibited, sets never intersect
            const auto expected_x_node =
                independent ? std::nullopt : std::optional< double >{ 0. };
            fractures.add_x_node_interaction(
                { "fset_01", "fset_02" }, beta, expected_x_node );

            run_and_validate( engine, fractures,
                absl::StrCat( "sim_two_fracture_sets_test_", beta ) );
        }
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_attractive_x_node_is_rejected()
    {
        geode::Logger::info( "TEST - X-NODE INTERACTION WITH BETA > 1" );

        // NOLINTBEGIN(*-magic-numbers)
        geode::FractureProcessBuilder fractures;
        auto& fset_01 = fractures.add_fracture_set( "fset_01", 0.05 );
        geode_unused( fset_01 );
        auto& fset_02 = fractures.add_fracture_set( "fset_02", 0.05 );
        geode_unused( fset_02 );
        bool rejected{ false };
        try
        {
            fractures.add_x_node_interaction( { "fset_01", "fset_02" }, 2. );
        }
        catch( const geode::OpenGeodeException& )
        {
            rejected = true;
        }
        // NOLINTEND(*-magic-numbers)
        geode::OpenGeodeStochasticStochasticException::test( rejected,
            "[Fractures] attractive X-node interaction (beta > 1) should be "
            "rejected" );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_fracture_density();
        test_fracture_intensity();
        test_fracture_minimal_spacing();
        test_observed_fractures();
        test_two_fracture_sets();
        test_attractive_x_node_is_rejected();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
