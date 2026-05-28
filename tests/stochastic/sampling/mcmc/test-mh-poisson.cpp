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
#include <geode/geometry/point.hpp>

#include <geode/stochastic/inference/statistics_tools.hpp>
#include <geode/stochastic/inference/target_statistics.hpp>

#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>
#include <geode/stochastic/models/gibbs_energy.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>
namespace
{
    struct SetDescription
    {
        std::string name;
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    using PoissonDensityDescription = geode::SingleObjectTermConfig;

    class PoissonConfig
    {
    public:
        PoissonConfig() = default;

        void add_domain_config( const geode::Point2D& min_p,
            const geode::Point2D& max_p,
            double buffer_size )
        {
            domain_config_.min_point = min_p;
            domain_config_.max_point = max_p;
            domain_config_.buffer_size = buffer_size;
        }

        void add_set_descriptor( const SetDescription& descriptor )
        {
            set_descriptors_.push_back( descriptor );
        }

        void add_density_descriptor(
            const PoissonDensityDescription& descriptor )
        {
            density_descriptors_.push_back( descriptor );
        }

        void add_target_statistics(
            const geode::TargetStatisticConfig& statistic_descriptor )
        {
            targeted_statistics_descriptors_.push_back( statistic_descriptor );
        }

        geode::SimulationContext< geode::Point2D > build() const
        {
            geode::SimulationContext< geode::Point2D > context;

            context.domain = geode::build_spatial_domain( domain_config_ );
            auto proposal_kernel = create_sets_and_set_samplers( context );
            create_model( context );

            context.mh_sampler =
                std::make_unique< geode::MetropolisHastings< geode::Point2D > >(
                    *context.model, std::move( proposal_kernel ) );
            create_target_statistics( context );
            return context;
        }

    private:
        std::unique_ptr< geode::ProposalKernel< geode::Point2D > >
            create_sets_and_set_samplers(
                geode::SimulationContext< geode::Point2D >& context ) const
        {
            auto proposal_kernel =
                std::make_unique< geode::ProposalKernel< geode::Point2D > >();
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id =
                    context.object_sets->add_set( set_desc.name );
                context.set_samplers.push_back(
                    std::make_unique< geode::UniformPointSetSampler< 2 > >(
                        *context.domain ) );
                geode::add_birth_death_change_moves(
                    context.set_samplers.back(), *proposal_kernel, set_id,
                    set_desc.birth_ratio, set_desc.death_ratio,
                    set_desc.change_ratio );
            }
            return proposal_kernel;
        }

        void create_model(
            geode::SimulationContext< geode::Point2D >& context ) const
        {
            geode::ModelConfig config;
            for( const auto& energy_desc : density_descriptors_ )
            {
                config.terms.emplace_back( energy_desc );
            }

            context.model = std::move( geode::build_model< geode::Point2D >(
                config, *context.object_sets, *context.domain ) );
        }

        void create_target_statistics(
            geode::SimulationContext< geode::Point2D >& context ) const
        {
            context.target_statistics.emplace( *context.model );
            for( const auto& target_stat : targeted_statistics_descriptors_ )
            {
                context.target_statistics->set_target( target_stat );
            }
        }

    private:
        geode::SpatialDomainConfig< 2 > domain_config_;
        std::vector< SetDescription > set_descriptors_;
        std::vector< PoissonDensityDescription > density_descriptors_;
        std::vector< geode::TargetStatisticConfig >
            targeted_statistics_descriptors_;
    };

    void test_single_type_poisson()
    {
        geode::Logger::info( "TEST - MH SINGLE TYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-POISSON@" );
        std::array< double, 4 > birth_ratio{ 0.1, 0.5, 2., 4. };
        std::array< double, 4 > change_ratio{ 0., 1., 1., 0. };

        for( const auto config : geode::Range{ birth_ratio.size() } )
        {
            PoissonConfig poisson_config;

            // NOLINTBEGIN(*-magic-numbers)
            poisson_config.add_domain_config( geode::Point2D{ { 0.0, 0.0 } },
                geode::Point2D{ { 10.0, 10.0 } }, 0. );
            // --- Set description
            SetDescription set_a;
            set_a.name = "A";
            set_a.birth_ratio = birth_ratio[config];
            set_a.death_ratio = 1.0;
            set_a.change_ratio = change_ratio[config];

            // --- Energy term description
            PoissonDensityDescription density_a;
            density_a.term_name = "density";
            density_a.object_set_names = { "A" };
            density_a.lambda = 0.3;
            density_a.object_feature = geode::ObjectInDomainFeatureConfig{};

            geode::TargetStatisticConfig stat_a{ "density", 30.0, 0.15 };

            poisson_config.add_set_descriptor( set_a );
            poisson_config.add_density_descriptor( density_a );
            poisson_config.add_target_statistics( stat_a );

            auto context = poisson_config.build();

            geode::SimulationRunner< geode::Point2D > runner{ std::move(
                context ) };

            // run simulation
            geode::SimulationPrinterConfigurator printer_config;
            printer_config.output_folder =
                absl::StrCat( printer_config.output_folder,
                    "/sim_point_poisson_test_", config );

            geode::SimulationConfigurator sim_config;
            sim_config.realizations = 1000;
            sim_config.metropolis_hasting_steps = 1000;
            sim_config.burn_in_steps = 1000;
            sim_config.printer = printer_config;

            auto statistic_tracker = runner.run( engine, sim_config );
            geode::statistics::validate(
                statistic_tracker, runner.target_statistics() );
        }
        // NOLINTEND(*-magic-numbers)

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_multitype_poisson()
    {
        geode::Logger::info( "TEST - MH MULTITYPE POISSON" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-POISSON-multi@" );
        PoissonConfig poisson_config;

        // NOLINTBEGIN(*-magic-numbers)
        poisson_config.add_domain_config( geode::Point2D{ { 0.0, 0.0 } },
            geode::Point2D{ { 10.0, 10.0 } }, 0. );

        // --- Set descriptions
        SetDescription set01{ "set01", 2.0, 3.0, 1.0 };
        SetDescription set02{ "set02", 3.0, 0.5, 1.0 };
        SetDescription set03{ "set03", 4.0, 1.0, 1.0 };

        // --- Energy term descriptions
        PoissonDensityDescription density01;
        density01.term_name = "density01";
        density01.object_set_names = { "set01" };
        density01.lambda = 0.1;
        density01.object_feature = geode::ObjectInDomainFeatureConfig{};

        geode::TargetStatisticConfig stat01{ "density01", 10.0, 0.15 };

        PoissonDensityDescription density02;
        density02.term_name = "density02";
        density02.object_set_names = { "set02" };
        density02.lambda = 0.4;
        density02.object_feature = geode::ObjectInDomainFeatureConfig{};

        geode::TargetStatisticConfig stat02{ "density02", 40.0, 0.15 };

        PoissonDensityDescription density03;
        density03.term_name = "density03";
        density03.object_set_names = { "set03" };
        density03.lambda = 0.3;
        density03.object_feature = geode::ObjectInDomainFeatureConfig{};

        geode::TargetStatisticConfig stat03{ "density03", 30.0, 0.15 };

        poisson_config.add_set_descriptor( set01 );
        poisson_config.add_set_descriptor( set02 );
        poisson_config.add_set_descriptor( set03 );

        poisson_config.add_density_descriptor( density01 );
        poisson_config.add_density_descriptor( density02 );
        poisson_config.add_density_descriptor( density03 );

        poisson_config.add_target_statistics( stat01 );
        poisson_config.add_target_statistics( stat02 );
        poisson_config.add_target_statistics( stat03 );

        geode::SimulationRunner< geode::Point2D > runner(
            poisson_config.build() );
        // run simulation
        geode::SimulationPrinterConfigurator printer_config;
        printer_config.output_folder = absl::StrCat(
            printer_config.output_folder, "/sim_point_multitype_poisson_test" );

        geode::SimulationConfigurator sim_config;
        sim_config.realizations = 1500;
        sim_config.metropolis_hasting_steps = 1000;
        sim_config.burn_in_steps = 3000;
        sim_config.printer = printer_config;
        // NOLINTEND(*-magic-numbers)

        auto statistic_tracker = runner.run( engine, sim_config );
        geode::statistics::validate(
            statistic_tracker, runner.target_statistics() );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_single_type_poisson();
        test_multitype_poisson();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}