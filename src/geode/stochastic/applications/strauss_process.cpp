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

#include <geode/stochastic/applications/strauss_process.hpp>

namespace
{
    std::vector< std::pair< std::string, std::string > >
        intra_set_interaction_names(
            const std::vector< std::string >& set_names )
    {
        std::vector< std::pair< std::string, std::string > > interaction_names;

        for( const auto name_id : geode::Range{ set_names.size() } )
        {
            interaction_names.emplace_back(
                set_names[name_id], set_names[name_id] );
        }
        return interaction_names;
    }

    std::vector< std::pair< std::string, std::string > >
        inter_set_interaction_names(
            const std::vector< std::string >& set_names )
    {
        std::vector< std::pair< std::string, std::string > > interaction_names;
        for( const auto id1 : geode::Range{ set_names.size() } )
        {
            for( const auto id2 : geode::Range{ id1 + 1, set_names.size() } )
            {
                interaction_names.emplace_back(
                    set_names[id1], set_names[id2] );
            }
        }
        return interaction_names;
    }

    std::string interaction_name(
        const std::vector< std::string >& set_names, double gamma )
    {
        std::string name{ "pwint_" };
        for( const auto& set_name : set_names )
        {
            absl::StrAppend( &name, set_name, "_" );
        }
        absl::StrAppend( &name, gamma );
        return name;
    }
} // namespace

namespace geode
{
    template < typename ObjectType >
    void StraussProcessBuilder< ObjectType >::set_domain(
        const SpatialDomainConfig< ObjectType::dim >& domain_cfg )
    {
        context_cfg_.domain = domain_cfg;
    };

    template < typename ObjectType >
    ObjectSetDefinition< ObjectType >&
        StraussProcessBuilder< ObjectType >::add_set( std::string_view name,
            double lambda,
            std::optional< double > expected_count )
    {
        auto& object_set_cfg = context_cfg_.add_set( name );

        SingleObjectTermConfig density;

        density.term_name = absl::StrCat( name, "_density" );
        density.object_set_names = { std::string( name ) };
        density.lambda = lambda;
        density.object_feature = ObjectInDomainFeatureConfig{};
        if( expected_count )
        {
            expected_stats_.push_back( geode::TargetStatisticConfig{
                density.term_name, *expected_count } );
        }
        context_cfg_.model.terms.emplace_back( std::move( density ) );

        return object_set_cfg;
    };

    template < typename ObjectType >
    void StraussProcessBuilder< ObjectType >::add_interaction(
        const std::vector< std::string >& interacting_set_names,
        double gamma,
        double distance_threshold,
        std::optional< double > expected_count,
        bool intra_set_interaction )
    {
        PairwiseTermConfig interaction;
        interaction.term_name =
            interaction_name( interacting_set_names, gamma );
        if( intra_set_interaction )
        {
            interaction.object_set_names_interactions =
                intra_set_interaction_names( interacting_set_names );
        }
        else
        {
            interaction.object_set_names_interactions =
                inter_set_interaction_names( interacting_set_names );
        }
        interaction.gamma = gamma;
        interaction.interaction_config =
            geode::MinimalDistanceCutoffConfig{ distance_threshold };

        if( expected_count )
        {
            expected_stats_.push_back( geode::TargetStatisticConfig{
                interaction.term_name, *expected_count } );
        }

        context_cfg_.model.terms.emplace_back( std::move( interaction ) );
    }

    template < typename ObjectType >
    SimulationContext< ObjectType >
        StraussProcessBuilder< ObjectType >::build_simulation_context() const
    {
        auto context =
            geode::build_simulation_context< ObjectType >( context_cfg_ );
        return context;
    };

    template < typename ObjectType >
    const std::vector< TargetStatisticConfig >&
        StraussProcessBuilder< ObjectType >::expected_statistics() const
    {
        return expected_stats_;
    };
    template class opengeode_stochastic_stochastic_api
        StraussProcessBuilder< Point2D >;
    template class opengeode_stochastic_stochastic_api
        StraussProcessBuilder< Point3D >;
    template class opengeode_stochastic_stochastic_api
        StraussProcessBuilder< geode::OwnerSegment2D >;

    //    using DensityDescription = geode::SingleObjectTermConfig;
    //    using InteractionDescription = geode::PairwiseTermConfig;
    //
    //    template < typename ObjectType >
    //    SimulationContext< ObjectType > build_strauss_process(
    //        const StraussProcessDescription< ObjectType >& desc )
    //    {
    //        SimulationContextConfig< ObjectType > config;
    //
    //        config.domain = desc.domain;
    //
    //        for( const auto& set_desc : desc.sets )
    //        {
    //            auto& set = config.add_set( set_desc.set_name );
    //
    //            set.sampler = set_desc.sampler;
    //
    //            set.dynamics.birth_ratio = set_desc.birth_ratio;
    //            set.dynamics.death_ratio = set_desc.death_ratio;
    //            set.dynamics.change_ratio = set_desc.change_ratio;
    //
    //            DensityDescription density;
    //
    //            density.term_name = set_desc.density_name;
    //            density.object_set_names = { set_desc.set_name };
    //            density.lambda = set_desc.lambda;
    //            density.object_feature = ObjectInDomainFeatureConfig{};
    //            config.model.terms.emplace_back( std::move( density ) );
    //        }
    //
    //        for( const auto& interaction_desc : desc.interactions )
    //        {
    //            InteractionDescription interaction;
    //            interaction.term_name = interaction_desc.interaction_name;
    //            interaction.object_set_names_interactions =
    //                std::move( build_interaction_set_names( interaction_desc )
    //                );
    //            interaction.gamma = interaction_desc.gamma;
    //            interaction.interaction_config =
    //                geode::MinimalDistanceCutoffConfig{
    //                interaction_desc.distance };
    //            config.model.terms.emplace_back( std::move( interaction ) );
    //        }
    //
    //        return build_simulation_context( config );
    //    }
    //
    //    template opengeode_stochastic_stochastic_api SimulationContext<
    //    Point2D >
    //        build_strauss_process< Point2D >(
    //            const StraussProcessDescription< Point2D >& );
    //    template opengeode_stochastic_stochastic_api SimulationContext<
    //    Point3D >
    //        build_strauss_process< Point3D >(
    //            const StraussProcessDescription< Point3D >& );
    //
    //    template < typename ObjectType >
    //    std::vector< geode::TargetStatisticConfig >
    //    build_strauss_targeted_stat(
    //        const StraussProcessDescription< ObjectType >& description )
    //    {
    //        std::vector< geode::TargetStatisticConfig > targets;
    //
    //        for( const auto& set_desc : description.sets )
    //        {
    //            if( !set_desc.expected_nb_objects )
    //            {
    //                continue;
    //            }
    //            targets.push_back( geode::TargetStatisticConfig{
    //                set_desc.density_name, *set_desc.expected_nb_objects } );
    //        }
    //        for( const auto& inter_desc : description.interactions )
    //        {
    //            if( !inter_desc.expected_nb_interactions )
    //            {
    //                continue;
    //            }
    //            targets.push_back(
    //                geode::TargetStatisticConfig{ inter_desc.interaction_name,
    //                    *inter_desc.expected_nb_interactions } );
    //        }
    //
    //        return targets;
    //    }
    //
    //    template opengeode_stochastic_stochastic_api
    //        std::vector< geode::TargetStatisticConfig >
    //        build_strauss_targeted_stat< Point2D >(
    //            const StraussProcessDescription< Point2D >& );
    //    template opengeode_stochastic_stochastic_api
    //        std::vector< geode::TargetStatisticConfig >
    //        build_strauss_targeted_stat< Point3D >(
    //            const StraussProcessDescription< Point3D >& );
} // namespace geode