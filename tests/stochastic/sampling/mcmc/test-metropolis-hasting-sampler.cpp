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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/models/energy_term_collection.hpp>
#include <geode/stochastic/models/energy_terms/density_term.hpp>
#include <geode/stochastic/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
namespace
{
    void test_acceptance_prob_helper()
    {
        // NOLINTBEGIN(*-magic-numbers)

        // log_accept >= 0 → prob = 1
        geode::OpenGeodeStochasticStochasticException::test(
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                0.5 )
                == 1.0,
            "[MH test] acceptance_prob_helper wrong for positive log_accept." );

        // very negative → prob = 0
        geode::OpenGeodeStochasticStochasticException::test(
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                -800.0 )
                == 0.0,
            "[MH test] acceptance_prob_helper wrong for extreme negative." );

        // moderate negative → exp(log_accept)
        auto val =
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                -1.0 );
        geode::OpenGeodeStochasticStochasticException::test(
            std::fabs( val - std::exp( -1.0 ) ) < geode::GLOBAL_EPSILON,
            "[MH test] acceptance_prob_helper wrong for -1.0." );
        // NOLINTEND(*-magic-numbers)
    }

    void test_beta_setter( geode::MetropolisHastings< geode::Point2D >& mh_eng )
    {
        // NOLINTBEGIN(*-magic-numbers)
        mh_eng.set_beta( 0.5 );
        geode::OpenGeodeStochasticStochasticException::test(
            mh_eng.beta() == 0.5, "[MH test] beta not set correctly." );

        bool exception_thrown = false;
        try
        {
            mh_eng.set_beta( -1.0 );
        }
        catch( ... )
        {
            exception_thrown = true;
        }
        geode::OpenGeodeStochasticStochasticException::test(
            exception_thrown, "[MH test] negative beta did not throw." );
        // NOLINTEND(*-magic-numbers)
    }

    void test_steps( const geode::MetropolisHastings< geode::Point2D >& mh_eng,
        geode::ObjectSets< geode::Point2D >& state )
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::RandomEngine engine;

        geode::index_t stat_sum{ 0 };
        constexpr geode::index_t NUNBER_ITR{ 100000 };

        geode::index_t accepted_birth{ 0 };
        geode::index_t accepted_death{ 0 };
        geode::index_t accepted_change{ 0 };
        geode::index_t nb_accepted{ 0 };

        for( const auto count : geode::Range{ NUNBER_ITR } )
        {
            auto result = mh_eng.step( state, engine );
            // Invariant: fixed object must remain

            geode::OpenGeodeStochasticStochasticException::test(
                result.decision == geode::MH_DECISION::accepted
                    || result.decision == geode::MH_DECISION::rejected,
                "[MH test] decision should be accepted or rejected." );

            // Log each step (optional: comment out if too verbose)
            //            geode::Logger::info( "Step: ", count,
            //                " move_type= ", static_cast< int >(
            //                result.move_type ), " decision= ", result.decision
            //                == geode::MH_DECISION::accepted ? "accepted"
            //                                                               :
            //                                                               "rejected",
            //                " delta_log_energy = ", result.delta_log_energy,
            //                " log_accept = ", result.log_accept,
            //                " state_size = ", state.size() );

            // Keep track of accepted moves by type
            if( result.decision == geode::MH_DECISION::accepted )
            {
                nb_accepted++;
                switch( result.move_type )
                {
                    case geode::MoveType::Birth:
                        accepted_birth++;
                        break;
                    case geode::MoveType::Death:
                        accepted_death++;
                        break;
                    case geode::MoveType::Change:
                        accepted_change++;
                        break;
                    default:
                        break;
                }
            }
            // should be change... only pone group here
            stat_sum += state.nb_objects();

            if( count % 1000 == 0 )
            {
                geode::Logger::info( "Progress: ", count, "  ", N,
                    " Mean objects =  ",
                    static_cast< double >( stat_sum )
                        / static_cast< double >( count ),
                    " nb accepted = ", nb_accepted, " accepted(B/D/C) = ",
                    static_cast< double >( accepted_birth )
                        / static_cast< double >( nb_accepted ),
                    "  ",
                    static_cast< double >( accepted_death )
                        / static_cast< double >( nb_accepted ),
                    "  ",
                    static_cast< double >( accepted_change )
                        / static_cast< double >( nb_accepted ) );
            }
        }
        // NOLINTEND(*-magic-numbers)
    }

} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        // NOLINTBEGIN(*-magic-numbers)
        geode::Point2D min_point{ { 0., 0. } };
        geode::Point2D max_point{ { 10., 10. } };

        geode::BoundingBox2D box;
        box.add_point( min_point );
        box.add_point( max_point );
        geode::SpatialDomain domain( box, 0. );

        geode::ObjectSets< geode::Point2D > state;
        const auto set_id = state.add_set( "default_name" );
        geode::UniformPointSetSampler< 2 > sampler( domain );

        double birth_prob = 0.3;
        double death_prob = 0.1;
        auto kernel = geode::create_birth_death_change_kernel< geode::Point2D >(
            set_id, sampler, birth_prob, death_prob );

        geode::EnergyTermCollection< geode::Point2D > energy_terms;

        // Add intensity term
        const auto term_id = energy_terms.add_energy_term(
            std::make_unique< geode::DensityTerm< geode::Point2D > >(
                "intensity", 0.5, std::vector< geode::uuid >{ set_id },
                domain ) );
        geode_unused( term_id );

        geode::MetropolisHastings< geode::Point2D > mh_eng(
            energy_terms, std::move( kernel ) );

        state.add_object( geode::Point2D{ { 1., 1. } }, set_id, true );
        state.add_object( geode::Point2D{ { 2., 2. } }, set_id, false );
        state.add_object( geode::Point2D{ { 3., 3. } }, set_id, true );

        test_steps( mh_eng, state );
        // geode::OpenGeodeStochasticStochasticException::test( state.get_set(
        // set_id ).nb_fixed_objects() == 2
        // );

        test_beta_setter( mh_eng );
        test_acceptance_prob_helper();

        geode::Logger::info( "MH TEST SUCCESS" );
        return 0;
        // NOLINTEND(*-magic-numbers)
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}