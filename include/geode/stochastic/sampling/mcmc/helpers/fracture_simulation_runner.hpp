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
#pragma once

#include <geode/basic/types.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_runner.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/density_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/spatial/pairwise_interactions.hpp>

namespace geode
{
    struct FractureSetDescription
    {
        std::string name;

        DoubleSampler::DistributionDescription length;
        DoubleSampler::DistributionDescription azimuth;

        // positionning
        double p20;
        // double p21;
        double minimal_spacing{ 0. };

        // mh dynamique
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };

        std::string string() const
        {
            auto message = absl::StrCat( "FractureSetDescription: ", name );
            absl::StrAppend(
                &message, "\n\t --> length distribution: ", length.string() );
            absl::StrAppend(
                &message, "\n\t --> azimuth distribution: ", azimuth.string() );
            absl::StrAppend( &message, "\n\t --> targeted p20: ", p20 );
            absl::StrAppend(
                &message, "\n\t --> minimal spacing: ", minimal_spacing );
            absl::StrAppend( &message,
                "\n\t --> MH move ratio - birth/death/change (", birth_ratio,
                " / ", death_ratio, " / ", change_ratio, ")" );
            return message;
        }
    };

    class FractureSimulationRunner : public SimulationRunner< OwnerSegment2D >
    {
    public:
        FractureSimulationRunner( const BoundingBox2D& box ) : box_( box ) {}

        void add_fracture_set_descriptor(
            const FractureSetDescription& descriptor )
        {
            set_descriptors_.push_back( descriptor );
        }

        void initialize() override
        {
            auto proposal_kernel =
                std::make_unique< ProposalKernel< OwnerSegment2D > >();

            // Mapping set names -> UUID
            std::unordered_map< std::string, uuid > name_to_uuid;

            // Step 1: create object sets and samplers
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = this->object_sets_.add_set( set_desc.name );
                name_to_uuid[set_desc.name] = set_id;

                auto length_distribution =
                    DoubleSampler::create_distribution( set_desc.length );
                auto azimuth_distribution =
                    DoubleSampler::create_angle_distribution_in_rad(
                        set_desc.azimuth );
                this->set_samplers_.push_back(
                    std::make_unique< UniformSegmentSetSampler >(
                        box_, length_distribution, azimuth_distribution ) );

                add_birth_death_change_moves( this->set_samplers_.back(),
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
                        std::make_unique< DensityTerm< OwnerSegment2D > >(
                            absl::StrCat( set_desc.name, "_density" ),
                            set_desc.p20,
                            absl::flat_hash_set< uuid >{ set_id } ) ) );
                // spacing
                if( set_desc.minimal_spacing < GLOBAL_EPSILON )
                {
                    continue;
                }
                auto interaction = std::make_unique<
                    EuclideanCutoffInteraction< OwnerSegment2D > >(
                    set_desc.minimal_spacing,
                    PairwiseInteraction< OwnerSegment2D >::SCOPE::same_set );

                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique< PairwiseTerm< OwnerSegment2D > >(
                            absl::StrCat( set_desc.name, "_min_spacing" ), 0.,
                            absl::flat_hash_set< uuid >{ set_id },
                            std::move( interaction ) ) ) );
            }

            this->mh_sampler_ =
                std::make_unique< MetropolisHastings< OwnerSegment2D > >(
                    this->energy_terms_collection_,
                    std::move( proposal_kernel ) );
        }

        void check_statistics(
            const StatisticsMonitor& statistic_monitoring ) const
        {
            const auto& computed_means = statistic_monitoring.means();

            for( const auto stat_id :
                Range{ this->energy_terms_collection_.size() } )
            {
                const auto& term = energy_terms_collection_.get(
                    ordered_energy_terms_[stat_id] );
                Logger::info( "[MH test] Statistic value ",
                    computed_means[stat_id],
                    " for energy term: ", term.name().data() );
            }
        }

    private:
        BoundingBox2D box_;
        std::vector< FractureSetDescription > set_descriptors_;
    };

} // namespace geode