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
            absl::string_view set_name, absl::string_view density_name )
        {
            auto& set = sets.emplace_back();
            set.set_name = set_name;
            set.density_name = density_name;
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