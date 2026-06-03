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

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
namespace
{
    std::vector< geode::Fracture > build_observed_fractures(
        const std::vector< std::array< geode::Point2D, 2 > >&
            fracture_extremities )
    {
        std::vector< geode::Fracture > fractures;
        for( const auto& extremities : fracture_extremities )
        {
            fractures.emplace_back(
                geode::OwnerSegment2D{ extremities[0], extremities[1] } );
        }
        return fractures;
    }
} // namespace

namespace geode
{
    using FractureDensityDescription = SingleObjectTermConfig;
    using FractureIntensityDescription = SingleObjectTermConfig;
    using FractureSpacingDescription = PairwiseTermConfig;

    using FractureSimulationConfig = SimulationContextConfig< Fracture >;

    FractureSimulationContext build_fractures_simulation_context(
        const FractureNetworkDescription& fnet_desc )
    {
        FractureSimulationConfig simulation_config;
        simulation_config.domain = fnet_desc.domain;

        for( const auto& fset_desc : fnet_desc.fracture_sets )
        {
            auto& fset = simulation_config.add_set( fset_desc.fset_name );

            fset.sampler = fset_desc.sampler;
            fset.dynamics.birth_ratio = fset_desc.birth_ratio;
            fset.dynamics.death_ratio = fset_desc.death_ratio;
            fset.dynamics.change_ratio = fset_desc.change_ratio;
            fset.fixed_objects =
                build_observed_fractures( fset_desc.observed_fractures );

            FractureDensityDescription density;
            density.term_name = fset_desc.density_name();
            density.object_set_names = { fset_desc.fset_name };
            density.lambda = fset_desc.p20;
            density.object_feature = ObjectInDomainFeatureConfig{};
            simulation_config.model.terms.emplace_back( std::move( density ) );

            FractureIntensityDescription intensity;
            intensity.term_name = fset_desc.intensity_name();
            intensity.object_set_names = { fset_desc.fset_name };
            intensity.lambda = fset_desc.p21;
            constexpr double caracteristic_length = 1.0;
            intensity.object_feature =
                SegmentLengthInsideBoxFeatureConfig{ caracteristic_length };
            simulation_config.model.terms.emplace_back(
                std::move( intensity ) );

            FractureSpacingDescription spacing;
            spacing.term_name = fset_desc.spacing_name();
            spacing.object_set_names_interactions = { { fset_desc.fset_name,
                fset_desc.fset_name } };
            spacing.gamma = 0.;
            spacing.interaction_config =
                geode::MinimalDistanceCutoffConfig{ fset_desc.minimal_spacing };
            simulation_config.model.terms.emplace_back( std::move( spacing ) );
        }

        return build_simulation_context( simulation_config );
    }

    std::vector< geode::TargetStatisticConfig > build_fractures_targeted_stat(
        const FractureNetworkDescription& description )
    {
        std::vector< geode::TargetStatisticConfig > targets;

        for( const auto& set_desc : description.fracture_sets )
        {
            if( set_desc.expected_number )
            {
                targets.push_back( geode::TargetStatisticConfig{
                    set_desc.density_name(), *set_desc.expected_number } );
            }
        }

        return targets;
    }

} // namespace geode