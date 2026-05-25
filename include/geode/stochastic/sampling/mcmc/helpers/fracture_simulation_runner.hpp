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
#pragma once

#include <geode/basic/types.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/models/energy_terms/density_term.hpp>
#include <geode/stochastic/models/energy_terms/intensity_term.hpp>
#include <geode/stochastic/models/energy_terms/pairwise_term.hpp>
#include <geode/stochastic/sampling/direct/double_sampler.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

#include <geode/geometry/basic_objects/segment.hpp>

namespace geode
{
    struct FractureSetDescription
    {
        std::string name;

        DoubleSampler::DistributionDescription length;
        DoubleSampler::DistributionDescription azimuth;

        // positionning
        double p20;
        std::vector< std::array< geode::Point2D, 2 > > observed_fractures;
        double p21{ 1. };
        double minimal_spacing{ 0. };

        // mh dynamique
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };

        std::string string() const
        {
            auto message = absl::StrCat( "FractureSetDescription: ", name );
            for( const auto& fixed_object : observed_fractures )
            {
                absl::StrAppend( &message,
                    "\n\t --> observation (x,y,z)start: ",
                    fixed_object[0].string(),
                    " (x,y,z)end: ", fixed_object[1].string() );
            }
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
        FractureSimulationRunner( const SpatialDomain< 2 >& domain )
            : SimulationRunner< OwnerSegment2D >( domain )
        {
        }

        void add_x_node_monitoring( double beta_x_node )
        {
            OpenGeodeStochasticStochasticException::check_exception(
                beta_x_node <= 1.0 && beta_x_node >= 0., nullptr,
                OpenGeodeException::TYPE::data,
                "[FractureSimulationRunner] x node should be inhibitated, "
                "please provise a value in [0., 1.]." );
            beta_x_node_ = beta_x_node;
        }
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
            //           std::unordered_map< std::string, uuid >
            //           set_name_to_uuid_;

            // Step 1: create object sets and samplers
            for( const auto& set_desc : set_descriptors_ )
            {
                const auto set_id = this->object_sets_.add_set( set_desc.name );
                for( const auto& fixed_object : set_desc.observed_fractures )
                {
                    this->object_sets_.add_object(
                        geode::OwnerSegment2D{
                            fixed_object[0], fixed_object[1] },
                        set_id, true );
                }
                set_name_to_uuid_[set_desc.name] = set_id;

                auto length_distribution =
                    DoubleSampler::create_distribution( set_desc.length );
                auto azimuth_distribution =
                    DoubleSampler::create_rad_angle_distribution_from_degree(
                        set_desc.azimuth );
                this->set_samplers_.push_back(
                    std::make_unique< UniformSegmentSetSampler >(
                        domain_, length_distribution, azimuth_distribution ) );

                add_birth_death_change_moves( this->set_samplers_.back(),
                    *proposal_kernel, set_id, set_desc.birth_ratio,
                    set_desc.death_ratio, set_desc.change_ratio );
            }

            // Step 2: create density energy terms
            for( const auto& set_desc : set_descriptors_ )
            {
                set_density_term( set_desc );
                set_intensity_term( set_desc );
                set_minimal_spacing_term( set_desc );
            }
            set_x_intersection_term();

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
                    computed_means[stat_id], " for energy term: ",
                    term.name().value_or( term.id().string() ) );
            }
        }

        std::string string() const
        {
            auto message =
                absl::StrCat( "Fracture Simulation Runner description" );
            for( const auto& desc : set_descriptors_ )
            {
                absl::StrAppend( &message, "\n\t ", desc.string() );
            }
            if( std::fabs( beta_x_node_ - 1. ) > GLOBAL_EPSILON )
            {
                absl::StrAppend( &message,
                    "\n\t --> x node monitioring (beta inhibition value): ",
                    beta_x_node_ );
            }
            else
            {
                absl::StrAppend(
                    &message, "\n\t --> x node monitioring : no inhibition." );
            }
            return message;
        }

    private:
        void set_density_term( const FractureSetDescription& set_desc )
        {
            const auto set_id = set_name_to_uuid_.at( set_desc.name );
            this->ordered_energy_terms_.push_back(
                this->energy_terms_collection_.add_energy_term(
                    std::make_unique< DensityTerm< OwnerSegment2D > >(
                        absl::StrCat( set_desc.name, "_density" ), set_desc.p20,
                        std::vector< uuid >{ set_id }, this->domain_ ) ) );
        }

        void set_intensity_term( const FractureSetDescription& set_desc )
        {
            if( std::fabs( set_desc.p21 - 1. ) < GLOBAL_EPSILON )
            {
                return;
            }
            const auto set_id = set_name_to_uuid_.at( set_desc.name );
            this->ordered_energy_terms_.push_back(
                this->energy_terms_collection_.add_energy_term(
                    std::make_unique< IntensityTerm >(
                        absl::StrCat( set_desc.name, "_intensity" ),
                        set_desc.p21, std::vector< uuid >{ set_id },
                        0.5 * this->domain_.smallest_length(),
                        this->domain_ ) ) );
        }

        void set_minimal_spacing_term( const FractureSetDescription& set_desc )
        {
            if( set_desc.minimal_spacing < GLOBAL_EPSILON )
            {
                return;
            }
            const auto set_id = set_name_to_uuid_.at( set_desc.name );
            auto interaction = std::make_unique<
                EuclideanCutoffInteraction< OwnerSegment2D > >(
                set_desc.minimal_spacing,
                PairwiseInteraction< OwnerSegment2D >::SCOPE::same_set );

            this->ordered_energy_terms_.push_back(
                this->energy_terms_collection_.add_energy_term(
                    std::make_unique< PairwiseTerm< OwnerSegment2D > >(
                        absl::StrCat( set_desc.name, "_min_spacing" ), 0.,
                        std::vector< uuid >{ set_id }, std::move( interaction ),
                        this->domain_ ) ) );
        }

        void set_x_intersection_term()
        {
            if( std::fabs( beta_x_node_ - 1. ) > GLOBAL_EPSILON )
            {
                std::vector< uuid > set_uuids;
                set_uuids.reserve( set_name_to_uuid_.size() );
                for( const auto& [name, id] : set_name_to_uuid_ )
                {
                    set_uuids.push_back( id );
                }
                auto interaction = std::make_unique<
                    EuclideanCutoffInteraction< OwnerSegment2D > >(
                    0., PairwiseInteraction<
                            OwnerSegment2D >::SCOPE::different_set );

                this->ordered_energy_terms_.push_back(
                    this->energy_terms_collection_.add_energy_term(
                        std::make_unique< PairwiseTerm< OwnerSegment2D > >(
                            absl::StrCat( "inter_set_x_nodes" ), beta_x_node_,
                            set_uuids, std::move( interaction ),
                            this->domain_ ) ) );
            }
        }

    private:
        std::vector< FractureSetDescription > set_descriptors_;
        std::unordered_map< std::string, uuid > set_name_to_uuid_;
        // x node monitoring
        double beta_x_node_{ 1. };
    };

} // namespace geode