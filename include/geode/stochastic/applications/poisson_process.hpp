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

#include <geode/stochastic/inference/target_statistics.hpp>

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>

namespace geode
{
    template < typename ObjectType >
    class PoissonProcessBuilder
    {
    public:
        void set_domain(
            const SpatialDomainConfig< ObjectType::dim >& domain_cfg );

        [[nodiscard]] ObjectSetDefinition< ObjectType >& add_set(
            std::string_view name,
            double lambda,
            std::optional< double > expected_count = std::nullopt );

        [[nodiscard]] SimulationContext< ObjectType >
            build_simulation_context() const;

        [[nodiscard]] const std::vector< TargetStatisticConfig >&
            expected_statistics() const;

    private:
        SimulationContextConfig< ObjectType > context_cfg_;
        std::vector< geode::TargetStatisticConfig > expected_stats_;
    };

} // namespace geode