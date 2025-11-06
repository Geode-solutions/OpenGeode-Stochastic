/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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
#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_runner.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/density_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

namespace
{
    struct FractureSetDescription
    {
        std::string name;

        geode::DoubleSampler::DistributionDescription length;
        geode::DoubleSampler::DistributionDescription azimuth;

        // positionning
        double p20;
        // double p21;
        double minimal_spacing;

        // mh dynamique
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    class FractureSimulationRunner
        : public geode::SimulationRunner< geode::OwnerSegment2D >
    {
    public:
        FractureSimulationRunner( const geode::BoundingBox2D& box )
            : box_( box )
        {
        }

        void add_fracture_set_descriptor(
            const FractureSetDescription& descriptor )
        {
            set_descriptors_.push_back( descriptor );
        }

        void initialize() override
        {
            auto proposal_kernel = std::make_unique<
                geode::ProposalKernel< geode::OwnerSegment2D > >();

            // Mapping set names -> UUID
            std::unordered_map< std::string, geode::uuid > name_to_uuid;

            // Step 1: create object sets and samplers
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = this->object_sets_.add_set( set_desc.name );
                name_to_uuid[set_desc.name] = set_id;

                auto length_distribution =
                    geode::DoubleSampler::create_distribution(
                        set_desc.length );
                auto azimuth_distribution =
                    geode::DoubleSampler::create_distribution(
                        set_desc.azimuth );
                this->set_samplers_.push_back(
                    std::make_unique< geode::UniformSegmentSetSampler >(
                        box_, length_distribution, azimuth_distribution ) );

                geode::add_birth_death_change_moves( this->set_samplers_.back(),
                    *proposal_kernel, set_id, set_desc.birth_ratio,
                    set_desc.death_ratio, set_desc.change_ratio );
            }

            // Step 2: create density energy terms
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = name_to_uuid.at( set_desc.name );
                // p20
                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique<
                            geode::DensityTerm< geode::OwnerSegment2D > >(
                            absl::StrCat( set_desc.name, "_density" ),
                            set_desc.p20,
                            absl::flat_hash_set< geode::uuid >{ set_id } ) ) );
                // spacing
                auto interaction =
                    std::make_unique< geode::EuclideanCutoffInteraction<
                        geode::OwnerSegment2D > >( set_desc.minimal_spacing,
                        geode::PairwiseInteraction<
                            geode::OwnerSegment2D >::SCOPE::same_set );

                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique<
                            geode::PairwiseTerm< geode::OwnerSegment2D > >(
                            absl::StrCat( set_desc.name, "_min_spacing" ), 0.,
                            absl::flat_hash_set< geode::uuid >{ set_id },
                            std::move( interaction ) ) ) );
            }

            this->mh_sampler_ = std::make_unique<
                geode::MetropolisHastings< geode::OwnerSegment2D > >(
                this->energy_terms_collection_, std::move( proposal_kernel ) );
        }

        void check_statistics(
            const geode::MonitoringStatistics& statistic_monitoring ) const
        {
            for( const auto stat_id :
                geode::Range{ this->energy_terms_collection_.size() } )
            {
                const auto& term = energy_terms_collection_.get(
                    ordered_energy_terms_[stat_id] );
                geode::Logger::info( "[MH test] Statistic value ",
                    statistic_monitoring.means[stat_id],
                    " for energy term: ", term.name().data() );
            }
        }

    private:
        geode::BoundingBox2D box_;
        std::vector< FractureSetDescription > set_descriptors_;
    };

    void test_fracture_simulator()
    {
        geode::Logger::info( "TEST - MH SINGLE SET FRACTURE SIMULATOR (with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-Fracture-set@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 100.0, 100.0 } } );

        // --- Object set
        FractureSetDescription setA;
        setA.name = "A";

        // length
        setA.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.length.min_value = 1;
        setA.length.max_value = 10.;

        // azimuth
        setA.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.azimuth.min_value = 1;
        setA.azimuth.max_value = 10.;

        // positionning
        setA.p20 = 0.1;
        setA.minimal_spacing = 1.;

        FractureSimulationRunner runner( box );
        runner.add_fracture_set_descriptor( setA );

        runner.initialize();

        constexpr geode::index_t steps = 1000;
        constexpr geode::index_t nb_realizations = 750;

        runner.run( engine, 10000 );
        auto stats = runner.run_print_and_monitor(
            "single_fracture_set_stats", engine, steps, nb_realizations );
        runner.check_statistics( stats );

        geode::Logger::info( "--> SUCCESS!" );
    }

    void test_two_fracture_sets_simulator()
    {
        geode::Logger::info( "TEST - MH TWO SET FRACTURE SIMULATOR (with "
                             "intra-set interactions)" );

        geode::RandomEngine engine;
        engine.set_seed( "@mh-test-single-Fracture-set@" );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 100.0, 100.0 } } );

        // --- Object set
        FractureSetDescription setA;
        setA.name = "A";

        // length
        setA.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.length.min_value = 1;
        setA.length.max_value = 10.;

        // azimuth
        setA.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setA.azimuth.min_value = 1;
        setA.azimuth.max_value = 10.;

        // positionning
        setA.p20 = 0.1;
        setA.minimal_spacing = 1.;

        // --- Object set
        FractureSetDescription setB;
        setB.name = "B";

        // length
        setB.length.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setB.length.min_value = 1;
        setB.length.max_value = 10.;

        // azimuth
        setB.azimuth.distribution_type =
            geode::UniformClosed< double >::distribution_type_static();
        setB.azimuth.min_value = 90.;
        setB.azimuth.max_value = 100.;

        // positionning
        setB.p20 = 0.1;
        setB.minimal_spacing = 2.;

        FractureSimulationRunner runner( box );
        runner.add_fracture_set_descriptor( setA );
        runner.add_fracture_set_descriptor( setB );

        runner.initialize();

        constexpr geode::index_t steps = 1000;
        constexpr geode::index_t nb_realizations = 750;

        runner.run( engine, 10000 );
        auto stats = runner.run_print_and_monitor(
            "two_fracture_sets_stats", engine, steps, nb_realizations );
        runner.check_statistics( stats );

        geode::Logger::info( "--> SUCCESS!" );
    }
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_fracture_simulator();
        test_two_fracture_sets_simulator();
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}