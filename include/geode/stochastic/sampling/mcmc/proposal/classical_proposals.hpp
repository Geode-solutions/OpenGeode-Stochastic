
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
#include <geode/stochastic/sampling/mcmc/proposal/marked_object_sampler/marked_object_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>

namespace geode
{
    template < typename Geometry >
    std::unique_ptr< ProposalKernel< Geometry > > create_birth_death_kernel(
        const MarkedObjectSampler< Geometry >& sampler,
        double birth_prob = 0.5,
        double death_prob = 0.5 )
    {
        auto kernel = std::make_unique< ProposalKernel< Geometry > >();
        kernel->add_move(
            std::make_unique< BirthMove< Geometry > >( sampler, birth_prob ) );
        kernel->add_move(
            std::make_unique< DeathMove< Geometry > >( sampler, death_prob ) );
        return kernel;
    }

    template < typename Geometry >
    std::unique_ptr< ProposalKernel< Geometry > >
        create_birth_death_change_kernel(
            const MarkedObjectSampler< Geometry >& sampler,
            double birth_prob = 0.33,
            double death_prob = 0.33,
            double change_prob = 0.34 )
    {
        auto kernel = std::make_unique< ProposalKernel< Geometry > >();
        kernel->add_move(
            std::make_unique< BirthMove< Geometry > >( sampler, birth_prob ) );
        kernel->add_move(
            std::make_unique< DeathMove< Geometry > >( sampler, death_prob ) );
        kernel->add_move( std::make_unique< ChangeMove< Geometry > >(
            sampler, change_prob ) );
        return kernel;
    }
} // namespace geode