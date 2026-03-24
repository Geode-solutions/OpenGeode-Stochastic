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

#include <geode/basic/assert.hpp>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/sampling/mcmc/energy_terms/density_term.hpp>
#include <geode/stochastic/sampling/mcmc/energy_terms/energy_term_collection.hpp>
#include <geode/stochastic/sampling/mcmc/energy_terms/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/energy_terms/pairwise_term.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

void test_gibbs_energy()
{
    geode::BoundingBox2D box;
    box.add_point( geode::Point2D{ { -1.0, -1.0 } } );
    box.add_point( geode::Point2D{ { 1.0, 1.0 } } );
    geode::SpatialDomain domain( box, 0.5 );

    geode::ObjectSets< geode::Point2D > pattern;
    const auto set_id = pattern.add_set( "default_name" );
    geode::Point2D p1{ { 0., 0. } };
    geode::Point2D p2{ { 1., 1. } };
    pattern.add_object( std::move( p1 ), set_id, false );
    pattern.add_object( std::move( p2 ), set_id, false );

    geode::EnergyTermCollection< geode::Point2D > energy_terms;

    // Add intensity term
    const auto term_id = energy_terms.add_energy_term(
        std::make_unique< geode::DensityTerm< geode::Point2D > >(
            "intensity", 0.5, std::vector< geode::uuid >{ set_id }, domain ) );
    geode_unused( term_id );
    // Add pairwise term with trivial interaction: always counts 1 for each pair
    auto interaction =
        std::make_unique< geode::EuclideanCutoffInteraction< geode::Point2D > >(
            1000000 );

    const auto pwterm_id = energy_terms.add_energy_term(
        std::make_unique< geode::PairwiseTerm< geode::Point2D > >(
            "interaction", 0.8, std::vector< geode::uuid >{ set_id },
            std::move( interaction ), domain ) );
    geode_unused( pwterm_id );

    OPENGEODE_EXCEPTION( energy_terms.size() == 2,
        "[test gibbs] Wrong number of components after adding terms." );

    geode::GibbsEnergy< geode::Point2D > gibbs_energy( energy_terms );

    // Check total log-energy is finite
    double total_energy = gibbs_energy.total_log_energy( pattern );
    OPENGEODE_EXCEPTION( std::isfinite( total_energy ),
        "[test gibbs] Total energy should be finite." );

    // Add new point to test delta_add
    geode::Point2D p3{ { 2., 2. } };
    geode::ObjectRef< geode::Point2D > p_ref{ p3, set_id };
    double delta_add = gibbs_energy.delta_log_add( pattern, p_ref );
    OPENGEODE_EXCEPTION( std::isfinite( delta_add ),
        "[test gibbs] Delta add should be finite." );

    geode::ObjectId obj_id{ 0, false, set_id };
    // Remove point test
    double delta_remove = gibbs_energy.delta_log_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION( std::isfinite( delta_remove ),
        "[test gibbs] Delta remove should be finite." );

    // Change point test
    double delta_change =
        gibbs_energy.delta_log_change( pattern, obj_id, p_ref );
    OPENGEODE_EXCEPTION( std::isfinite( delta_change ),
        "[test gibbs] Delta change should be finite." );

    // Clear components and verify
    energy_terms.clear();
    OPENGEODE_EXCEPTION( energy_terms.size() == 0,
        "[test gibbs] Components not cleared properly." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        test_gibbs_energy();
        geode::Logger::info( "MH TEST GIBBS ENERGY SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
    return 0;
}