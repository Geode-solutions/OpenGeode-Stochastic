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
    struct PoissonSetDescription
    {
        std::string set_name;

        ObjectSamplerConfig< ObjectType > sampler;

        std::string density_name;
        double lambda{ 0. };
        std::optional< double > expected_nb_objects;

        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    template < typename ObjectType >
    struct PoissonProcessDescription
    {
        std::string process_name{ "Poisson" };
        SpatialDomainConfig< ObjectType::dim > domain;

        std::vector< PoissonSetDescription< ObjectType > > sets;

        PoissonSetDescription< ObjectType >& add_set(
            absl::string_view set_name )
        {
            auto& set = sets.emplace_back();
            set.set_name = set_name;
            set.density_name = absl::StrCat( set_name, "_density" );
            return set;
        }
    };

    template < typename ObjectType >
    SimulationContext< ObjectType > build_poisson_process(
        const PoissonProcessDescription< ObjectType >& description );

    template < typename ObjectType >
    std::vector< geode::TargetStatisticConfig > build_poisson_targeted_stat(
        const PoissonProcessDescription< ObjectType >& description );

} // namespace geode