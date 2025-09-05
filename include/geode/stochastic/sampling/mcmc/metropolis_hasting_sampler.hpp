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

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>

namespace geode
{
    enum struct MHDecision
    {
        Accepted,
        Rejected,
        Undecided
    };

    template < typename Object >
    struct StepResult
    {
        MHDecision decision{ MHDecision::Undecided };
        typename Proposal< Object >::Type move_type{
            Proposal< Object >::Type::Invalid
        };
        double log_accept{ -std::numeric_limits< double >::infinity() };
        double delta_log_energy{ 0.0 };
    };

    template < typename Object >
    class MetropolisHastings
    {
    public:
        MetropolisHastings( GibbsEnergy< Object >& energy,
            std::unique_ptr< ProposalKernel< Object > > proposal_kernel )
            : energy_( energy ),
              proposal_kernel_( std::move( proposal_kernel ) )
        {
            OPENGEODE_ASSERT(
                proposal_kernel_ != nullptr, "[MH] null proposal kernel" );
        }

        Configuration< Object > initialise_configuration_with_sampling(
            RandomEngine& engine,
            const std::unordered_map< uuid, index_t >& group_targets ) const
        {
            Configuration< Object > config;
            for( const auto& [gid, target] : group_targets )
            {
                config.add_group( gid );
                while( config.nb_objects_in_group( gid ) < target )
                {
                    bool added = false;
                    for( int attempt = 0; attempt < 100 && !added; ++attempt )
                    {
                        Proposal< Object > proposal =
                            proposal_kernel_->propose( config, engine );

                        if( proposal.type == Proposal< Object >::Type::Birth
                            && proposal.new_object->second == gid )
                        {
                            config.add_object(
                                std::move( proposal.new_object->first ),
                                proposal.new_object->second );
                            added = true;
                        }
                    }
                    OPENGEODE_EXCEPTION( added,
                        "[MH] Birth move need to be more probable for group: ",
                        gid.string() );
                }
            }
            return config;
        }

        //       Configuration< Object > initialise_configuration_with_burnin(
        //           RandomEngine& engine, index_t number_of_steps ) const
        //       {
        //           Configuration< Object > configuration;
        //
        //           for( const auto count : geode::Range{ number_of_steps } )
        //           {
        //               geode_unused( count );
        //               step( configuration, engine );
        //           }
        //           return configuration;
        //       }

        StepResult< Object > step(
            Configuration< Object >& state, RandomEngine& engine ) const
        {
            Proposal< Object > proposal =
                proposal_kernel_->propose( state, engine );

            if( proposal.type == Proposal< Object >::Type::Birth )
            {
                return birth_step( proposal, state, engine );
            }
            if( proposal.type == Proposal< Object >::Type::Death )
            {
                return death_step( proposal, state, engine );
            }
            if( proposal.type == Proposal< Object >::Type::Change )
            {
                return change_step( proposal, state, engine );
            }
            return StepResult< Object >{};
        }

        void walk( Configuration< Object >& state,
            RandomEngine& engine,
            index_t nb_steps ) const
        {
            for( const auto count : geode::Range{ nb_steps } )
            {
                geode_unused( count );
                step( state, engine );
            }
        }

        Configuration< Object > walk_copy( Configuration< Object > initial,
            RandomEngine& engine,
            index_t nb_steps ) const
        {
            walk( initial, engine, nb_steps );
            return initial;
        }

        double beta() const
        {
            return beta_;
        }

        void set_beta( double b )
        {
            OPENGEODE_EXCEPTION( b >= 0.0, "[MH] beta must be >= 0" );
            if( b == 0 )
            {
                geode::Logger::info(
                    "[Metropolis Hastings] - beta == 0 all move will be "
                    "accepted - Uniform sampling." );
            }
            if( b < 1 )
            {
                geode::Logger::info(
                    "[Metropolis Hastings] - beta < 1 moves that increase "
                    "energy are "
                    "more likely to be accepted - Hot system introduce "
                    "randomness for exploration." );
            }
            if( b == 1 )
            {
                geode::Logger::info( "[Metropolis Hastings] - beta == 1 "
                                     "default choice no temperature "
                                     "- only consider energy." );
            }
            if( b > 1 )
            {
                geode::Logger::info(
                    "[Metropolis Hastings] - beta > 1 moves that increase "
                    "energy are "
                    "less likely to be accepted - Cold system to ensure "
                    "convergence but may find local minimum randomness." );
            }
            beta_ = b;
        }

        static double acceptance_prob_helper( double log_accept )
        {
            if( std::isnan( log_accept ) )
                return 0.0;
            if( log_accept >= 0.0 )
                return 1.0;
            // prevent expoential overflow
            constexpr double LOG_MIN = -745.0;
            if( log_accept < LOG_MIN )
                return 0.0;
            return std::exp( log_accept );
        }

    private:
        const double compute_log_accept(
            const double deltaU, const Proposal< Object >& proposal ) const
        {
            OPENGEODE_ASSERT(
                std::isfinite( proposal.log_forward_prob )
                    && std::isfinite( proposal.log_backward_prob ),
                "[MH] Non-finite proposal log-probabilities" );
            return -beta_ * deltaU + proposal.log_backward_prob
                   - proposal.log_forward_prob;
        }

        template < typename ApplyMove >
        StepResult< Object > accept_or_reject( Proposal< Object >& proposal,
            Configuration< Object >& state,
            RandomEngine& engine,
            const double delta_log_energy,
            ApplyMove&& apply_move ) const
        {
            StepResult< Object > step_result;
            step_result.move_type = proposal.type;
            step_result.delta_log_energy = delta_log_energy;
            step_result.log_accept =
                compute_log_accept( delta_log_energy, proposal );

            double log_u = engine.sample_log();
            step_result.decision = ( log_u < step_result.log_accept )
                                       ? MHDecision::Accepted
                                       : MHDecision::Rejected;
            if( step_result.decision == MHDecision::Accepted )
                apply_move( state, proposal );

            return step_result;
        }

        StepResult< Object > birth_step( Proposal< Object >& proposal,
            Configuration< Object >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.new_object.has_value(),
                "[MH] Birth proposal has no new_object" );
            const auto delta_log_energy = energy_.delta_log_energy_add( state,
                proposal.new_object.value().first,
                proposal.new_object.value().second );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.add_object( std::move( p.new_object.value().first ),
                        p.new_object.value().second );
                } );
        };

        StepResult< Object > death_step( Proposal< Object >& proposal,
            Configuration< Object >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.old_object_id.has_value(),
                "[MH] Death proposal has no index" );
            const auto delta_log_energy = energy_.delta_log_energy_remove(
                state, proposal.old_object_id.value() );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.remove_object( p.old_object_id.value() );
                } );
        };

        StepResult< Object > change_step( Proposal< Object >& proposal,
            Configuration< Object >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.new_object.has_value(),
                "[MH] Change proposal has no new_object" );
            OPENGEODE_ASSERT( proposal.old_object_id.has_value(),
                "[MH] Change proposal has no index" );
            const auto delta_log_energy = energy_.delta_log_energy_change(
                state, proposal.old_object_id.value(),
                proposal.new_object.value().first );
            // should we test that objects are in the same group?
            // should be ensured by the dynamic
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.update_object( p.old_object_id.value(),
                        std::move( p.new_object.value().first ) );
                } );
        };

    private:
        const GibbsEnergy< Object >& energy_;
        std::unique_ptr< ProposalKernel< Object > > proposal_kernel_;
        double beta_{ 1.0 };
    };
} // namespace geode