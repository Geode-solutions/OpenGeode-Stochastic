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
namespace geode
{
    template < typename Geometry >
    class MetropolisHastings
    {
    public:
        MetropolisHastings( GibbsEnergy< Geometry >& energy,
            std::unique_ptr< ProposalKernel< Geometry > > proposal_kernel,
            double beta = 1.0 )
            : energy_( energy ),
              proposal_kernel_( std::move( proposal_kernel ) ),
              beta_( beta )
        {
        }

        void step( Configuration< Geometry >& state, std::mt19937& rng )
        {
            // propose a move
            Proposal< Geometry > proposal = proposal_kernel_->propose( state );

            // compute ΔU (energy difference)
            double deltaU = 0.0;
            switch( proposal.type )
            {
                case Proposal< Geometry >::Type::Birth:
                    deltaU = energy_.log_energy_delta_add(
                        state, *proposal.new_object );
                    break;
                case Proposal< Geometry >::Type::Death:
                    deltaU = energy_.log_energy_delta_remove(
                        state, *proposal.index );
                    break;
                case Proposal< Geometry >::Type::Change:
                    deltaU = energy_.log_energy_delta_change(
                        state, *proposal.index, *proposal.new_object );
                    break;
            }

            // acceptance probability
            double log_accept = -beta_ * deltaU + proposal.log_backward_prob
                                - proposal.log_forward_prob;

            std::uniform_real_distribution< double > dist( 0.0, 1.0 );
            double u = std::log( dist( rng ) );

            if( u < log_accept )
            {
                // accept move
                switch( proposal.type )
                {
                    case Proposal< Geometry >::Type::Birth:
                        state.add_object( std::move( *proposal.new_object ) );
                        break;
                    case Proposal< Geometry >::Type::Death:
                        state.remove_object( *proposal.index );
                        break;
                    case Proposal< Geometry >::Type::Change:
                        state.change_object( *proposal.index,
                            std::move( *proposal.new_object ) );
                        break;
                }
            }
            else
            {
                // reject → do nothing
            }
        }

    private:
        GibbsEnergy< Geometry >& energy_;
        std::unique_ptr< ProposalKernel< Geometry > > proposal_kernel_;
        double beta_;
    };
} // namespace geode