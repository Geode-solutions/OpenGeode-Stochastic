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

#include <geode/stochastic/sampling/distributions.hpp>

#include <geode/basic/pimpl.hpp>

namespace geode
{
    class opengeode_stochastic_stochastic_api RandomEngine
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( RandomEngine );

    public:
        explicit RandomEngine();
        ~RandomEngine();

        void set_seed( uint64_t number );
        void set_seed( std::string_view word );

        template < typename Type >
        Type sample_uniform( const UniformClosed< Type > &law );

        template < typename Type >
        Type sample_uniform( const UniformClosedOpen< Type > &law );

        double sample_gaussian( const Gaussian &law );
        double sample_truncated_gaussian( const TruncatedGaussian &law );

        double sample_von_mises( const VonMises &law );

        double sample_truncated_lognormal( const TruncatedLogNormal &law );
        double sample_truncated_powerlaw( const TruncatedPowerLaw &law );

        double sample_log();

        bool sample_bernoulli( double probability_of_success );

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };
} // namespace geode