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

#include <geode/stochastic/applications/observed_intensity.hpp>
#include <geode/stochastic/inference/target_statistics.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace geode
{
    using Fracture = OwnerSegment2D;
    using FractureSamplerConfig = ObjectSamplerConfig< Fracture >;
    using FractureSimulationContext = SimulationContext< Fracture >;
    using FractureSimulationRunner = SimulationRunner< Fracture >;

    class opengeode_stochastic_stochastic_api FractureProcessBuilder
    {
    public:
        void set_domain( const SpatialDomainConfig< 2 >& domain_cfg );

        /// Fracture set with a density term: p20 is the number per unit area
        /// of the fractures whose length follows the length law of the set
        /// (sampler.length). Its statistic counts the fractures in the domain
        /// by their endpoints (half the number of endpoints strictly inside),
        /// whichever endpoint anchors a fracture.
        ///
        /// Observed fractures go to the fixed_objects of the returned
        /// definition. Set the length law, p20 and scope from what was
        /// observed in the domain:
        /// - a sample of the fractures of all lengths: length law and p20 of
        /// all the fractures, IntensityScope::all_objects;
        /// - all the fractures longer than a cutoff L_obs, none shorter:
        /// length law cut at L_obs (max_value = L_obs), p20 of the fractures
        /// shorter than L_obs, IntensityScope::all_objects;
        /// - all the fractures longer than L_obs and a sample of the shorter
        /// ones: same settings as above.
        /// With all_objects, the observed fractures within the length law are
        /// a part of p20, and only the missing ones are simulated; the observed
        /// fractures outside it (longer than L_obs) are added as they are.
        /// This is exact without interaction terms on the set, and an
        /// approximation with them. In all the cases, the observed fractures
        /// act on the simulated ones through the interaction terms of the
        /// set.
        /// With IntensityScope::free_objects, p20 is the density of the
        /// simulated fractures only, and all the observed ones are added.
        ///
        /// expected_count, if given, is the expected mean of the statistic:
        /// p20 * |D| with all_objects, plus the count of the observed
        /// fractures outside the length law (or of all of them with
        /// free_objects).
        [[nodiscard]] ObjectSetDefinition< Fracture >& add_fracture_set(
            std::string_view name,
            double p20,
            std::optional< double > expected_count = std::nullopt,
            IntensityScope scope = IntensityScope::free_objects );

        /// Intensity term (fracture length per unit area) on a fracture set.
        void add_intensity( std::string_view set_name,
            double p21,
            std::optional< double > expected_total_length = std::nullopt );

        /// Hard-core constraint: fractures of the set cannot be closer than
        /// the given distance (0 forbids intersections).
        void add_minimal_spacing( std::string_view set_name,
            double minimal_spacing,
            std::optional< double > expected_count = std::nullopt );

        /// Interaction between intersecting fractures of different sets
        /// (beta in [0, 1]: 0 forbids X-nodes, 1 has no effect).
        void add_x_node_interaction(
            const std::vector< std::string >& interacting_set_names,
            double beta,
            std::optional< double > expected_count = std::nullopt );

        [[nodiscard]] FractureSimulationContext
            build_simulation_context() const;

        [[nodiscard]] const std::vector< TargetStatisticConfig >&
            expected_statistics() const;

    private:
        SimulationContextConfig< Fracture > context_cfg_;
        std::vector< TargetStatisticConfig > expected_stats_;
        // indices in sets of the sets in IntensityScope::all_objects, with
        // their p20
        std::vector< std::pair< index_t, double > > all_objects_sets_;
    };
} // namespace geode