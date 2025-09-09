
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
#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>

namespace geode
{
    template < typename Type >
    std::unique_ptr< ProposalKernel< Type > > create_birth_death_kernel(
        const ObjectSetSampler< Type >& sampler, double birth_prob )
    {
        auto kernel = std::make_unique< ProposalKernel< Type > >();
        kernel->add_move( std::make_unique< BirthDeathMove< Type > >(
            sampler, 1., birth_prob ) );
        return kernel;
    }

    template < typename Type >
    std::unique_ptr< ProposalKernel< Type > > create_birth_death_change_kernel(
        const ObjectSetSampler< Type >& sampler,
        double birth_prob,
        double death_prob )
    {
        auto birth_death_prob = birth_prob + death_prob;
        OPENGEODE_EXCEPTION( birth_death_prob < 1.,
            "[Proposal Kernel] - changes should be allowed." );
        auto kernel = std::make_unique< ProposalKernel< Type > >();
        kernel->add_move( std::make_unique< BirthDeathMove< Type > >(
            sampler, birth_death_prob, birth_prob / birth_death_prob ) );
        kernel->add_move( std::make_unique< ChangeMove< Type > >(
            sampler, 1. - birth_death_prob ) );
        return kernel;
    }
} // namespace geode