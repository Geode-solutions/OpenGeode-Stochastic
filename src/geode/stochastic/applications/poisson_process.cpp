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
    template < typename ObjectType >
    void PoissonProcessBuilder< ObjectType >::set_domain(
        const SpatialDomainConfig< ObjectType::dim >& domain_cfg )
    {
        context_cfg_.domain = domain_cfg;
    };

    template < typename ObjectType >
    ObjectSetDefinition< ObjectType >&
        PoissonProcessBuilder< ObjectType >::add_set( std::string_view name,
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
    SimulationContext< ObjectType >
        PoissonProcessBuilder< ObjectType >::build_simulation_context() const
    {
        auto context =
            geode::build_simulation_context< ObjectType >( context_cfg_ );
        return context;
    };

    template < typename ObjectType >
    const std::vector< TargetStatisticConfig >&
        PoissonProcessBuilder< ObjectType >::expected_statistics() const
    {
        return expected_stats_;
    };
    template class opengeode_stochastic_stochastic_api
        PoissonProcessBuilder< Point2D >;
    template class opengeode_stochastic_stochastic_api
        PoissonProcessBuilder< Point3D >;
    template class opengeode_stochastic_stochastic_api
        PoissonProcessBuilder< geode::OwnerSegment2D >;

} // namespace geode