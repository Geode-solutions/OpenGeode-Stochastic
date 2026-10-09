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

#include <absl/strings/str_join.h>

#include <geode/stochastic/applications/fractures.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
namespace
{
    std::vector< std::pair< std::string, std::string > > inter_set_interactions(
        const std::vector< std::string >& set_names )
    {
        std::vector< std::pair< std::string, std::string > > interactions;
        for( const auto id1 : geode::Range{ set_names.size() } )
        {
            for( const auto id2 : geode::Range{ id1 + 1, set_names.size() } )
            {
                interactions.emplace_back( set_names[id1], set_names[id2] );
            }
        }
        return interactions;
    }
} // namespace

namespace geode
{
    using FractureDensityDescription = SingleObjectTermConfig;
    using FractureIntensityDescription = SingleObjectTermConfig;
    using FractureSpacingDescription = PairwiseTermConfig;

    using XNodeInteractionDescription = geode::PairwiseTermConfig;

    void FractureProcessBuilder::set_domain(
        const SpatialDomainConfig< 2 >& domain_cfg )
    {
        context_cfg_.domain = domain_cfg;
    }

    ObjectSetDefinition< Fracture >& FractureProcessBuilder::add_fracture_set(
        std::string_view name,
        double p20,
        std::optional< double > expected_count )
    {
        auto& fset_cfg = context_cfg_.add_set( name );

        FractureDensityDescription density;
        density.term_name = absl::StrCat( name, "_p20" );
        density.object_set_names = { std::string( name ) };
        density.lambda = p20;
        density.object_feature = ObjectInDomainFeatureConfig{};
        if( expected_count )
        {
            expected_stats_.push_back(
                TargetStatisticConfig{ density.term_name, *expected_count } );
        }
        context_cfg_.model.terms.emplace_back( std::move( density ) );

        return fset_cfg;
    }

    void FractureProcessBuilder::add_intensity( std::string_view set_name,
        double p21,
        std::optional< double > expected_total_length )
    {
        FractureIntensityDescription intensity;
        intensity.term_name = absl::StrCat( set_name, "_p21" );
        intensity.object_set_names = { std::string( set_name ) };
        intensity.lambda = p21;
        constexpr double CARACTERISTIC_LENGTH = 1.0;
        intensity.object_feature =
            SegmentLengthInsideBoxFeatureConfig{ CARACTERISTIC_LENGTH };
        if( expected_total_length )
        {
            expected_stats_.push_back( TargetStatisticConfig{
                intensity.term_name, *expected_total_length } );
        }
        context_cfg_.model.terms.emplace_back( std::move( intensity ) );
    }

    void FractureProcessBuilder::add_minimal_spacing( std::string_view set_name,
        double minimal_spacing,
        std::optional< double > expected_count )
    {
        FractureSpacingDescription spacing;
        spacing.term_name = absl::StrCat( set_name, "_spacing" );
        spacing.object_set_names_interactions = { { std::string( set_name ),
            std::string( set_name ) } };
        spacing.gamma = 0.;
        spacing.interaction_config =
            MinimalDistanceCutoffConfig{ minimal_spacing };
        if( expected_count )
        {
            expected_stats_.push_back(
                TargetStatisticConfig{ spacing.term_name, *expected_count } );
        }
        context_cfg_.model.terms.emplace_back( std::move( spacing ) );
    }

    void FractureProcessBuilder::add_x_node_interaction(
        const std::vector< std::string >& interacting_set_names,
        double beta,
        std::optional< double > expected_count )
    {
        OpenGeodeStochasticStochasticException::check_exception(
            beta <= 1.0 && beta >= 0., nullptr, OpenGeodeException::TYPE::data,
            "[FractureProcessBuilder] x node should be inhibited, please "
            "provide a value in [0., 1.]." );
        OpenGeodeStochasticStochasticException::check_exception(
            interacting_set_names.size() > 1, nullptr,
            OpenGeodeException::TYPE::data,
            "[FractureProcessBuilder] x node interaction requires at least two "
            "fracture sets." );

        XNodeInteractionDescription interaction;
        interaction.term_name = absl::StrCat(
            "x_node_", absl::StrJoin( interacting_set_names, "_" ) );
        interaction.object_set_names_interactions =
            inter_set_interactions( interacting_set_names );
        interaction.gamma = beta;
        interaction.interaction_config = MinimalDistanceCutoffConfig{ 0. };
        if( expected_count )
        {
            expected_stats_.push_back( TargetStatisticConfig{
                interaction.term_name, *expected_count } );
        }
        context_cfg_.model.terms.emplace_back( std::move( interaction ) );
    }

    FractureSimulationContext
        FractureProcessBuilder::build_simulation_context() const
    {
        return geode::build_simulation_context< Fracture >( context_cfg_ );
    }

    const std::vector< TargetStatisticConfig >&
        FractureProcessBuilder::expected_statistics() const
    {
        return expected_stats_;
    }

} // namespace geode