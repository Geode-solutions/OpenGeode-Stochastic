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
        UnDecided
    };

    template < typename Geometry >
    struct StepResult
    {
        MHDecision decision{ MHDecision::UnDecided };
        typename Proposal< Geometry >::Type move_type{
            Proposal< Geometry >::Type::Invalid
        };
        double log_accept{ -std::numeric_limits< double >::infinity() };
        double delta_log_energy{ 0.0 };
    };

    template < typename Geometry >
    class MetropolisHastings
    {
    public:
        MetropolisHastings( GibbsEnergy< Geometry >& energy,
            std::unique_ptr< ProposalKernel< Geometry > > proposal_kernel )
            : energy_( energy ),
              proposal_kernel_( std::move( proposal_kernel ) )
        {
            OPENGEODE_ASSERT(
                proposal_kernel_ != nullptr, "[MH] null proposal kernel" );
        }

        StepResult< Geometry > step(
            Configuration< Geometry >& state, RandomEngine& engine ) const
        {
            Proposal< Geometry > proposal =
                proposal_kernel_->propose( state, engine );

            if( proposal.type == Proposal< Geometry >::Type::Birth )
            {
                return birth_step( proposal, state, engine );
            }
            if( proposal.type == Proposal< Geometry >::Type::Death )
            {
                return death_step( proposal, state, engine );
            }
            if( proposal.type == Proposal< Geometry >::Type::Change )
            {
                return change_step( proposal, state, engine );
            }
            return StepResult< Geometry >{};
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
                geode::Logger::info( "[MH] - beta == 0 all move will be "
                                     "accepted - Uniform sampling." );
            }
            if( b < 1 )
            {
                geode::Logger::info(
                    "[MH] - beta < 1 moves that increase energy are "
                    "more likely to be accepted - Hot system introduce "
                    "randomness for exploration." );
            }
            if( b == 1 )
            {
                geode::Logger::info( "[MH] - beta == 1 Standars choice "
                                     "- only consider energy." );
            }
            if( b > 1 )
            {
                geode::Logger::info(
                    "[MH] - beta > 1 moves that increase energy are "
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
            const double deltaU, const Proposal< Geometry >& proposal ) const
        {
            OPENGEODE_ASSERT( std::isfinite( log_forward_prob )
                                  && std::isfinite( log_backward_prob ),
                "[MH] Non-finite proposal log-probabilities" );
            return -beta_ * deltaU + proposal.log_backward_prob
                   - proposal.log_forward_prob;
        }

        template < typename ApplyMove >
        StepResult< Geometry > accept_or_reject( Proposal< Geometry >& proposal,
            Configuration< Geometry >& state,
            RandomEngine& engine,
            const double delta_log_energy,
            ApplyMove&& apply_move ) const
        {
            StepResult< Geometry > step_result;
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

        StepResult< Geometry > birth_step( Proposal< Geometry >& proposal,
            Configuration< Geometry >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.new_object.has_value(),
                "[MH] Birth proposal has no new_object" );
            const auto delta_log_energy = energy_.delta_log_energy_add(
                state, proposal.new_object.value() );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.add_object( std::move( p.new_object.value() ) );
                } );
        };

        StepResult< Geometry > death_step( Proposal< Geometry >& proposal,
            Configuration< Geometry >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.index.has_value(),
                "[MH] Death proposal has no index" );
            const auto delta_log_energy = energy_.delta_log_energy_remove(
                state, proposal.index.value() );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.remove_object( p.index.value() );
                } );
        };

        StepResult< Geometry > change_step( Proposal< Geometry >& proposal,
            Configuration< Geometry >& state,
            RandomEngine& engine ) const
        {
            OPENGEODE_ASSERT( proposal.new_object.has_value(),
                "[MH] Change proposal has no new_object" );
            OPENGEODE_ASSERT( proposal.index.has_value(),
                "[MH] Change proposal has no index" );
            const auto delta_log_energy = energy_.delta_log_energy_change(
                state, proposal.index.value(), proposal.new_object.value() );
            return accept_or_reject( proposal, state, engine, delta_log_energy,
                []( auto& s, auto& p ) {
                    s.change_object(
                        p.index.value(), std::move( p.new_object.value() ) );
                } );
        };

    private:
        const GibbsEnergy< Geometry >& energy_;
        std::unique_ptr< ProposalKernel< Geometry > > proposal_kernel_;
        double beta_{ 1.0 };
    };
} // namespace geode