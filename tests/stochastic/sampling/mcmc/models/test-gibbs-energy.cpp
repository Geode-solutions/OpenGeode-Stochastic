/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/stochastic/sampling/mcmc/models/components/intensity_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/pairwise_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/spatial/object_set.hpp>

namespace
{
    geode::ObjectSet< geode::Point2D > create_object_set(
        const geode::uuid& subset_id )
    {
        geode::Point2D p1{ { 0., 0. } };
        geode::Point2D p2{ { 1., 1. } };

        geode::ObjectSet< geode::Point2D > pattern;
        pattern.add_object( std::move( p1 ), subset_id );
        pattern.add_object( std::move( p2 ), subset_id );

        return pattern;
    }
} // namespace

void test_gibbs_energy( const geode::uuid& subset_id )
{
    geode::GibbsEnergy< geode::Point2D > gibbs_energy;

    // Add intensity term
    gibbs_energy.add_energy_term(
        std::make_unique< geode::IntensityTerm< geode::Point2D > >(
            0.5, subset_id ) );

    // Add pairwise term with trivial interaction: always counts 1 for each pair
    auto interaction_fn = []( const geode::Point2D& a,
                              const geode::Point2D& b ) {
        geode_unused( a );
        geode_unused( b );
        return true;
    };
    gibbs_energy.add_energy_term( std::make_unique<
        geode::PairwiseTerm< geode::Point2D, decltype( interaction_fn ) > >(
        0.8, interaction_fn ) );

    OPENGEODE_EXCEPTION( gibbs_energy.number_of_energy_terms() == 2,
        "[test gibbs] Wrong number of components after adding terms." );

    auto pattern = create_object_set( subset_id );

    // Check total log-energy is finite
    double total_energy = gibbs_energy.total_log_energy( pattern );
    OPENGEODE_EXCEPTION( std::isfinite( total_energy ),
        "[test gibbs] Total energy should be finite." );

    // Add new point to test delta_add
    geode::Point2D p3{ { 2., 2. } };
    double delta_add =
        gibbs_energy.delta_log_energy_add( pattern, p3, subset_id );
    OPENGEODE_EXCEPTION( std::isfinite( delta_add ),
        "[test gibbs] Delta add should be finite." );

    geode::ObjectId obj_id{ 0, subset_id };
    // Remove point test
    double delta_remove =
        gibbs_energy.delta_log_energy_remove( pattern, obj_id );
    OPENGEODE_EXCEPTION( std::isfinite( delta_remove ),
        "[test gibbs] Delta remove should be finite." );

    // Change point test
    double delta_change =
        gibbs_energy.delta_log_energy_change( pattern, obj_id, p3 );
    OPENGEODE_EXCEPTION( std::isfinite( delta_change ),
        "[test gibbs] Delta change should be finite." );

    // Check string outputs don’t crash
    auto comp_str = gibbs_energy.ordered_energy_term_parameter_string();
    auto stats_str =
        gibbs_energy.ordered_energy_term_statistics_string( pattern );
    OPENGEODE_EXCEPTION(
        !comp_str.empty(), "[test gibbs] Components string is empty." );
    OPENGEODE_EXCEPTION(
        !stats_str.empty(), "[test gibbs] Statistics string is empty." );

    // Clear components and verify
    gibbs_energy.clear_energy_terms();
    OPENGEODE_EXCEPTION( gibbs_energy.number_of_energy_terms() == 0,
        "[test gibbs] Components not cleared properly." );
}

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::uuid subset_id;

        // test_gibbs_energy( subset_id );
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
    return 0;
}