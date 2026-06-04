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

#include <geode/stochastic/applications/poisson_process.hpp>
#include <geode/stochastic/inference/target_statistics.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>

namespace geode
{
    template < typename ObjectType >
    struct StraussInteractionDescription
    {
        std::string interaction_name;

        std::vector< std::string > set_names;

        double gamma{ 1. };
        double distance{ 0. };

        bool include_intra_set{ true };
        bool include_inter_set{ false };

        std::optional< double > expected_nb_interactions;
    };

    template < typename ObjectType >
    struct StraussProcessDescription
    {
        SpatialDomainConfig< ObjectType::dim > domain;

        std::vector< PoissonSetDescription< ObjectType > > sets;

        std::vector< StraussInteractionDescription< ObjectType > > interactions;

        PoissonSetDescription< ObjectType >& add_set(
            absl::string_view set_name )
        {
            auto& set = sets.emplace_back();
            set.set_name = set_name;
            set.density_name = absl::StrCat( set_name, "_density" );

            return set;
        }

        StraussInteractionDescription< ObjectType >& add_interaction(
            absl::string_view interaction_name )
        {
            auto& interaction = interactions.emplace_back();
            interaction.interaction_name = interaction_name;
            return interaction;
        }
    };

    template < typename ObjectType >
    SimulationContext< ObjectType > build_strauss_process(
        const StraussProcessDescription< ObjectType >& description );

    template < typename ObjectType >
    std::vector< geode::TargetStatisticConfig > build_strauss_targeted_stat(
        const StraussProcessDescription< ObjectType >& description );

} // namespace geode