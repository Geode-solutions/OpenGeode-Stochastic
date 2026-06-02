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

#include <geode/stochastic/applications/poisson_process.hpp>

namespace geode
{
    using PoissonDensityDescription = geode::SingleObjectTermConfig;

    template < typename ObjectType >
    SimulationContext< ObjectType > build_poisson_process(
        const PoissonProcessDescription< ObjectType >& desc )
    {
        SimulationContextConfig< ObjectType > config;

        config.domain = desc.domain;

        for( const auto& set_desc : desc.sets )
        {
            auto& set = config.add_set( set_desc.set_name );

            set.sampler = set_desc.sampler;

            set.dynamics.birth_ratio = set_desc.birth_ratio;
            set.dynamics.death_ratio = set_desc.death_ratio;
            set.dynamics.change_ratio = set_desc.change_ratio;

            PoissonDensityDescription density;

            density.term_name = set_desc.density_name;
            density.object_set_names = { set_desc.set_name };
            density.lambda = set_desc.lambda;
            density.object_feature = ObjectInDomainFeatureConfig{};

            config.model.terms.emplace_back( std::move( density ) );
        }

        return build_simulation_context( config );
    }

    template opengeode_stochastic_stochastic_api SimulationContext< Point2D >
        build_poisson_process< Point2D >(
            const PoissonProcessDescription< Point2D >& );
    template opengeode_stochastic_stochastic_api SimulationContext< Point3D >
        build_poisson_process< Point3D >(
            const PoissonProcessDescription< Point3D >& );

    template < typename ObjectType >
    std::vector< geode::TargetStatisticConfig > build_poisson_targeted_stat(
        const PoissonProcessDescription< ObjectType >& description )
    {
        std::vector< geode::TargetStatisticConfig > targets;

        for( const auto& set_desc : description.sets )
        {
            if( !set_desc.expected_nb_objects )
            {
                continue;
            }
            targets.push_back( geode::TargetStatisticConfig{
                set_desc.density_name, *set_desc.expected_nb_objects } );
        }

        return targets;
    }

    template opengeode_stochastic_stochastic_api
        std::vector< geode::TargetStatisticConfig >
        build_poisson_targeted_stat< Point2D >(
            const PoissonProcessDescription< Point2D >& );
    template opengeode_stochastic_stochastic_api
        std::vector< geode::TargetStatisticConfig >
        build_poisson_targeted_stat< Point3D >(
            const PoissonProcessDescription< Point3D >& );
} // namespace geode