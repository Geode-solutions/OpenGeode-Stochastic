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
#include <geode/stochastic/sampling/mcmc/proposal/moves.hpp>

namespace geode
{

    template < typename ObjectType >
    class ProposalKernel
    {
    public:
        virtual ~ProposalKernel() = default;

        Proposal< ObjectType > propose(
            const ObjectSet< ObjectType >& current, RandomEngine& engine ) const
        {
            OPENGEODE_EXCEPTION( !moves_.empty(),
                "[MCMC Proposal Kernel] - no move are defined in the Kernel." );
            if( cumulative_probs_.size() == 1 )
            {
                return moves_[0]->propose_move( current, engine );
            }
            auto rnd = engine.sample_uniform( uniform_closed_double_ );
            for( const auto proba_id : Range{ cumulative_probs_.size() } )
            {
                if( rnd <= cumulative_probs_[proba_id] )
                {
                    return moves_[proba_id]->propose_move( current, engine );
                }
            }
            OPENGEODE_ASSERT_NOT_REACHED(
                "[MCMC Proposal Kernel]: Should not be reached move pdf is "
                "correctly set." );
            return moves_.back()->propose_move( current, engine );
        }

        void add_move( std::unique_ptr< Move< ObjectType > > move )
        {
            moves_.push_back( std::move( move ) );
            compute_cumulative_sum_probs();
        }

    private:
        void compute_cumulative_sum_probs()
        {
            const auto n = moves_.size();
            cumulative_probs_.resize( n );
            if( n == 1 )
            {
                cumulative_probs_[0] = 1.;
                return;
            }
            double sum{ 0. };
            for( const auto prob_id : Range{ n } )
            {
                sum += moves_[prob_id]->probability();
                cumulative_probs_[prob_id] = sum;
            }
            auto total = cumulative_probs_.back();
            OPENGEODE_EXCEPTION( total > GLOBAL_EPSILON,
                "[MCMC Proposal Kernel] - Total "
                "probability is zero in Kernel." );
            std::transform( cumulative_probs_.begin(), cumulative_probs_.end(),
                cumulative_probs_.begin(), [total]( double p ) {
                    return p / total;
                } );
            cumulative_probs_.back() = 1.0;
        }

    private:
        std::vector< std::unique_ptr< Move< ObjectType > > > moves_;
        std::vector< double > cumulative_probs_;

        geode::UniformClosed< double > uniform_closed_double_;
    };
} // namespace geode