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
#include <geode/stochastic/sampling/mcmc/proposal/moves.hpp>

namespace geode
{
    template < typename ObjectType >
    struct Proposal
    {
        uuid set_id;

        MoveResult< ObjectType > proposed_move;
        ObjectRef< ObjectType > new_object()
        {
            OpenGeodeStochasticStochasticException::check_exception(
                proposed_move.new_object.has_value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[Proposal] Proposal has no new_object" );
            return ObjectRef< ObjectType >{ proposed_move.new_object.value(),
                set_id };
        };

        ObjectId old_object_id()
        {
            OpenGeodeStochasticStochasticException::check_exception(
                proposed_move.old_object_id.has_value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[Proposal] Proposal has no old_object_id" );
            return ObjectId{ proposed_move.old_object_id.value(), false,
                set_id };
        };

        std::string string() const
        {
            return absl::StrCat( "Move proposal on subset: ", set_id, " -- ",
                proposed_move.string() );
        }
    };

    template < typename ObjectType >
    class ProposalKernel
    {
    public:
        virtual ~ProposalKernel() = default;

        Proposal< ObjectType > propose( const ObjectSets< ObjectType >& current,
            RandomEngine& engine ) const
        {
            OpenGeodeStochasticStochasticException::check_exception(
                !set_moves_.empty(), nullptr, OpenGeodeException::TYPE::data,
                "[MCMC Proposal Kernel] - no move are defined in the Kernel." );
            auto rnd = engine.sample_uniform( uniform_distribution_closed_ );
            for( const auto proba_id : Range{ cumulative_probs_.size() } )
            {
                if( rnd <= cumulative_probs_[proba_id] )
                {
                    auto& [set_id, move] = set_moves_[proba_id];
                    return Proposal< ObjectType >{ set_id,
                        move->propose_move(
                            current.get_set( set_id ), engine ) };
                }
            }
            throw OpenGeodeStochasticStochasticException{ nullptr,
                OpenGeodeException::TYPE::internal,
                "[MCMC Proposal Kernel]: Should not be reached move pdf is "
                "correctly set." };
        }

        void add_move(
            const uuid& set_id, std::unique_ptr< Move< ObjectType > > move )
        {
            set_moves_.push_back( { set_id, std::move( move ) } );
            initialize_probabilities();
        }

        std::string string() const
        {
            auto message = absl::StrCat( "Proposal Kernel:",
                "\n\t - number of moves: ", set_moves_.size() );
            absl::StrAppend(
                &message, "\n\t --> move cumulative probabilities:" );
            for( const auto cumsum : cumulative_probs_ )
            {
                absl::StrAppend( &message, " ", cumsum );
            }
            for( const auto& [set_id, move] : set_moves_ )
            {
                absl::StrAppend( &message, " \n\t --> move on subset ",
                    set_id.string(), ": ", move->string() );
            }
            return message;
        }

    private:
        std::vector< double > compute_probabilities() const
        {
            std::vector< double > probabilities( set_moves_.size(), 0. );

            // Extract weights
            std::transform( set_moves_.begin(), set_moves_.end(),
                probabilities.begin(), []( const auto& move ) {
                    return move.second->proportion_weight();
                } );

            // Compute total
            const double total = std::accumulate(
                probabilities.begin(), probabilities.end(), 0.0 );

            // Ensure total > 0
            OpenGeodeStochasticStochasticException::check_exception(
                total > GLOBAL_EPSILON, nullptr,
                OpenGeodeException::TYPE::internal,
                "[MCMC Proposal Kernel] - Total "
                "probability is zero in Kernel." );

            // Normalize
            std::transform( probabilities.begin(), probabilities.end(),
                probabilities.begin(), [total]( double p ) {
                    return p / total;
                } );
            return probabilities;
        }

        void compute_cumulative_probabilities(
            const std::vector< double >& probabilities )
        {
            cumulative_probs_.resize( probabilities.size() );

            // Compute cumulative sum; works safely for empty or single-element
            // vectors
            std::partial_sum( probabilities.begin(), probabilities.end(),
                cumulative_probs_.begin() );

            if( !cumulative_probs_.empty() )
                cumulative_probs_.back() = 1.0; // ensure exact 1.0
        }
        void initialize_move_probabilities(
            const std::vector< double >& probabilities )
        {
            for( const auto move_id : geode::Range{ probabilities.size() } )
            {
                set_moves_[move_id].second->initialize_probability(
                    probabilities[move_id] );
            }
        }

        void initialize_probabilities()
        {
            auto probabilities = compute_probabilities();
            compute_cumulative_probabilities( probabilities );
            initialize_move_probabilities( probabilities );
        }

    private:
        std::vector< std::pair< uuid, std::unique_ptr< Move< ObjectType > > > >
            set_moves_;
        std::vector< double > cumulative_probs_;

        geode::UniformClosed< double > uniform_distribution_closed_;
    };
} // namespace geode