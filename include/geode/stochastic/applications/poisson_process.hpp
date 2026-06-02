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
        double lambda;
        std::optional< double > expected_nb_objects;

        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };
    };

    template < typename ObjectType >
    struct PoissonProcessDescription
    {
        SpatialDomainConfig< ObjectType::dim > domain;

        std::vector< PoissonSetDescription< ObjectType > > sets;

        PoissonSetDescription< ObjectType >& add_set(
            absl::string_view set_name, absl::string_view density_name )
        {
            auto& set = sets.emplace_back();
            set.set_name = set_name;
            set.density_name = density_name;
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