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

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/models/model.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>

namespace geode
{
    enum struct MHDecision
    {
        Accepted,
        Rejected,
        Undecided
    };

    template < typename ObjectType >
    struct StepResult
    {
        MHDecision decision{ MHDecision::Undecided };
        MoveType move_type{ MoveType::Invalid };
        double log_accept{ -std::numeric_limits< double >::infinity() };
        double delta_log_energy{ 0.0 };
    };

    template < typename ObjectType >
    class MetropolisHastings
    {
    public:
        MetropolisHastings( const Model< ObjectType >& model,
            std::unique_ptr< ProposalKernel< ObjectType > > proposal_kernel )
            : model_{ model }, proposal_kernel_( std::move( proposal_kernel ) )
        {
            OpenGeodeStochasticStochasticException::check_exception(
                proposal_kernel_ != nullptr, nullptr,
                OpenGeodeException::TYPE::data,
                "[MetropolisHastings] Proposal kernel is not defined." );
        }

        StepResult< ObjectType > step(
            ObjectSets< ObjectType >& state, RandomEngine& engine ) const
        {
            Proposal< ObjectType > proposal =
                proposal_kernel_->propose( state, engine );
            const auto& move_type = proposal.proposed_move.type;
            if( move_type == MoveType::Birth )
            {
                return birth_step( proposal, state, engine );
            }
            if( move_type == MoveType::Death )
            {
                return death_step( proposal, state, engine );
            }
            if( move_type == MoveType::Change )
            {
                return change_step( proposal, state, engine );
            }
            return StepResult< ObjectType >{};
        }

        void walk( ObjectSets< ObjectType >& state,
            RandomEngine& engine,
            index_t nb_steps ) const
        {
            for( const auto count : geode::Range{ nb_steps } )
            {
                geode_unused( count );
                step( state, engine );
            }
        }

        ObjectSets< ObjectType > walk_copy( ObjectSets< ObjectType > initial,
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
            OpenGeodeStochasticStochasticException::check_exception( b >= 0.0,
                nullptr, OpenGeodeException::TYPE::data,
                "[MetropolisHastings] The teperature (beta) must be >= 0" );
            if( b == 0 )
            {
                Logger::info( "[MetropolisHastings] beta == 0 all move will be "
                              "accepted - Uniform sampling." );
            }
            if( b < 1 )
            {
                Logger::info(
                    "[MetropolisHastings] beta < 1 moves that increase "
                    "energy are more likely to be accepted - Hot system "
                    "introduce randomness for exploration." );
            }
            if( b == 1 )
            {
                Logger::info( "[MetropolisHastings] beta == 1 default "
                              "choice no temperature - only consider energy." );
            }
            if( b > 1 )
            {
                Logger::info( "[MetropolisHastings] beta > 1 moves that "
                              "increase energy are less likely to be accepted "
                              "- Cold system to ensure convergence but may "
                              "find local minimum randomness." );
            }
            beta_ = b;
        }

        static double acceptance_prob_helper( double log_accept )
        {
            if( std::isnan( log_accept ) )
                return 0.0;
            if( log_accept >= 0.0 )
                return 1.0;
            // prevent exponential overflow
            constexpr double LOG_MIN = -745.0;
            if( log_accept < LOG_MIN )
                return 0.0;
            return std::exp( log_accept );
        }

    private:
        const double compute_log_accept( const double deltaU,
            const ProposalProbabilities& proposal_probas ) const
        {
            return -beta_ * deltaU + proposal_probas.transition_probability();
        }

        template < typename ApplyMove >
        StepResult< ObjectType > accept_or_reject(
            Proposal< ObjectType >& proposal,
            ObjectSets< ObjectType >& state,
            RandomEngine& engine,
            const double delta_log_energy,
            ApplyMove&& apply_move ) const
        {
            const auto& proposed_move = proposal.proposed_move;
            StepResult< ObjectType > step_result;
            step_result.move_type = proposed_move.type;
            step_result.delta_log_energy = delta_log_energy;
            step_result.log_accept = compute_log_accept(
                delta_log_energy, proposed_move.proposal_probabilities );

            double log_u = engine.sample_log();
            step_result.decision = ( log_u < step_result.log_accept )
                                       ? MHDecision::Accepted
                                       : MHDecision::Rejected;
            if( step_result.decision == MHDecision::Accepted )
                apply_move( state, proposal );

            return step_result;
        }

        StepResult< ObjectType > birth_step( Proposal< ObjectType >& proposal,
            ObjectSets< ObjectType >& state,
            RandomEngine& engine ) const
        {
            const auto new_object = proposal.new_object();
            const auto delta_log_energy =
                model_.energy().delta_log_add( state, new_object );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& state, auto& proposal ) {
                    state.add_object(
                        std::move( proposal.proposed_move.new_object.value() ),
                        proposal.set_id, false );
                } );
        };

        StepResult< ObjectType > death_step( Proposal< ObjectType >& proposal,
            ObjectSets< ObjectType >& state,
            RandomEngine& engine ) const
        {
            const auto old_object_id = proposal.old_object_id();
            const auto delta_log_energy =
                model_.energy().delta_log_remove( state, old_object_id );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& state, auto& proposal ) {
                    state.remove_free_object( proposal.old_object_id() );
                } );
        };

        StepResult< ObjectType > change_step( Proposal< ObjectType >& proposal,
            ObjectSets< ObjectType >& state,
            RandomEngine& engine ) const
        {
            const auto new_object = proposal.new_object();
            const auto old_object_id = proposal.old_object_id();
            const auto delta_log_energy = model_.energy().delta_log_change(
                state, old_object_id, new_object );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& state, auto& proposal ) {
                    state.update_free_object( proposal.old_object_id(),
                        std::move(
                            proposal.proposed_move.new_object.value() ) );
                } );
        };

    private:
        const Model< ObjectType >& model_;
        std::unique_ptr< ProposalKernel< ObjectType > > proposal_kernel_;
        double beta_{ 1.0 };
    };
} // namespace geode