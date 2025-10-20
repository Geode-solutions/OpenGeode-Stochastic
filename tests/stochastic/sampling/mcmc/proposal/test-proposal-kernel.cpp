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
#include <memory>

#include <geode/basic/assert.hpp>
#include <geode/geometry/point.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/moves.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace
{
    geode::uuid init_object_set( geode::ObjectSets< geode::Point2D >& pattern )
    {
        geode::Point2D p1{ { 0., 0. } };
        geode::Point2D p2{ { 1., 1. } };

        auto set_id = pattern.add_set();
        pattern.add_object( std::move( p1 ), set_id );
        pattern.add_object( std::move( p2 ), set_id );

        return set_id;
    }

    void test_proposal_kernel()
    {
        geode::ObjectSets< geode::Point2D > config;
        auto set_id = init_object_set( config );

        geode::Point2D min_point{ { 0., 0. } };
        geode::Point2D max_point{ { 10., 100. } };

        geode::BoundingBox2D box;
        box.add_point( min_point );
        box.add_point( max_point );

        geode::UniformPointSetSampler< 2 > sampler( box );

        // Create classical birth-death-change kernel
        auto kernel = geode::create_birth_death_change_kernel< geode::Point2D >(
            set_id, sampler, 0.4, 0.4 );

        geode::RandomEngine engine;

        bool saw_birth = false, saw_death = false, saw_change = false;

        for( const auto i : geode::Range{ 400 } )
        {
            auto proposal = kernel->propose( config, engine );
            const auto& proposed_move = proposal.proposed_move;

            switch( proposed_move.type )
            {
                case geode::MoveType::Birth:
                    saw_birth = true;
                    OPENGEODE_EXCEPTION( proposed_move.new_object.has_value(),
                        "[test proposal] Birth must provide new_object." );
                    OPENGEODE_EXCEPTION(
                        !proposed_move.old_object_id.has_value(),
                        "[test proposal] Birth should not provide index." );
                    // Probabilities
                    OPENGEODE_EXCEPTION(
                        proposed_move.proposal_probabilities.log_forward_prob
                            <= 0.0,
                        "[test proposal] Birth forward log-prob must be <= "
                        "0." );

                    OPENGEODE_EXCEPTION(
                        std::abs(
                            proposed_move.proposal_probabilities
                                .log_backward_prob
                            - ( std::log( 0.8 * 0.5 )
                                - std::log( config.nb_objects_in_set( set_id )
                                            + 1.0 ) ) )
                            < geode::GLOBAL_EPSILON,
                        "[test proposal] Birth backward log-prob mismatch." );
                    break;

                case geode::MoveType::Death:
                    saw_death = true;
                    OPENGEODE_EXCEPTION( !proposed_move.new_object.has_value(),
                        "[test proposal] Death should not provide "
                        "new_object." );
                    OPENGEODE_EXCEPTION(
                        proposed_move.old_object_id.has_value(),
                        "[test proposal] Death must provide index." );
                    OPENGEODE_EXCEPTION(
                        proposed_move.old_object_id.value()
                            < config.nb_objects_in_set( set_id ),
                        "[test proposal] Death index out of bounds." );
                    // Probabilities
                    OPENGEODE_EXCEPTION(
                        proposed_move.proposal_probabilities.log_forward_prob
                            <= 0.0,
                        "[test proposal] Death forward log-prob must be <= "
                        "0." );
                    break;

                case geode::MoveType::Change:
                    saw_change = true;
                    OPENGEODE_EXCEPTION( proposed_move.new_object.has_value(),
                        "[test proposal] Change must provide new_object." );
                    OPENGEODE_EXCEPTION(
                        proposed_move.old_object_id.has_value(),
                        "[test proposal] Change must provide index." );
                    OPENGEODE_EXCEPTION(
                        proposed_move.old_object_id.value()
                            < config.nb_objects_in_set( set_id ),
                        "[test proposal] Change index out of bounds." );
                    break;

                case geode::MoveType::Invalid:
                default:
                    throw geode::OpenGeodeException(
                        "[test proposal] Proposal type Invalid." );
            }

            // Log probabilities must be finite proposed move should be
            // possible.
            OPENGEODE_EXCEPTION(
                std::isfinite(
                    proposed_move.proposal_probabilities.log_forward_prob ),
                "[test proposal] Forward probability is not finite - Move" );
            OPENGEODE_EXCEPTION(
                std::isfinite(
                    proposed_move.proposal_probabilities.log_backward_prob ),
                "[test proposal] Backward probability is not finite." );
        }

        // Ensure kernel actually produced all types of proposals
        OPENGEODE_EXCEPTION( saw_birth && saw_death && saw_change,
            "[test proposal] Kernel did not produce all move types." );

        // --- Edge case: empty object_set ---
        geode::ObjectSets< geode::Point2D > empty_config;
        const auto empty_set_id = empty_config.add_set();
        auto empty_kernel =
            geode::create_birth_death_change_kernel< geode::Point2D >(
                empty_set_id, sampler, 0.4, 0.4 );
        auto proposal = empty_kernel->propose( empty_config, engine );
        const auto& proposed_move = proposal.proposed_move;

        OPENGEODE_EXCEPTION(
            proposed_move.type == geode::MoveType::Birth
                || proposed_move.type == geode::MoveType::Invalid,
            "[test proposal] On empty config, only Birth should be possible." );
    }
} // namespace
int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        test_proposal_kernel();
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
    return 0;
}