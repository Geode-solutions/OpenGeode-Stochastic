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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/sampling/direct/configuration_sampler/point_configuration_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/intensity_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
namespace
{
    void test_acceptance_prob_helper()
    {
        // log_accept >= 0 → prob = 1
        OPENGEODE_EXCEPTION(
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                0.5 )
                == 1.0,
            "[MH test] acceptance_prob_helper wrong for positive log_accept." );

        // very negative → prob = 0
        OPENGEODE_EXCEPTION(
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                -800.0 )
                == 0.0,
            "[MH test] acceptance_prob_helper wrong for extreme negative." );

        // moderate negative → exp(log_accept)
        double val =
            geode::MetropolisHastings< geode::Point2D >::acceptance_prob_helper(
                -1.0 );
        OPENGEODE_EXCEPTION( std::abs( val - std::exp( -1.0 ) ) < 1e-12,
            "[MH test] acceptance_prob_helper wrong for -1.0." );
    }

    void test_beta_setter( geode::MetropolisHastings< geode::Point2D >& mh )
    {
        mh.set_beta( 0.5 );
        OPENGEODE_EXCEPTION(
            mh.beta() == 0.5, "[MH test] beta not set correctly." );

        bool exception_thrown = false;
        try
        {
            mh.set_beta( -1.0 );
        }
        catch( ... )
        {
            exception_thrown = true;
        }
        OPENGEODE_EXCEPTION(
            exception_thrown, "[MH test] negative beta did not throw." );
    }

    void test_steps( const geode::MetropolisHastings< geode::Point2D >& mh,
        const geode::uuid& group_id )
    {
        geode::RandomEngine engine;

        std::unordered_map< geode::uuid, geode::index_t > targets = {
            { group_id, 20 }
        };
        geode::Configuration< geode::Point2D > state =
            mh.initialize_configuration_with_sampling( engine, targets );

        geode::index_t stat_sum{ 0 };
        constexpr geode::index_t N{ 100000 };

        geode::index_t accepted_birth{ 0 };
        geode::index_t accepted_death{ 0 };
        geode::index_t accepted_change{ 0 };
        geode::index_t nb_accepted{ 0 };

        for( const auto count : geode::Range{ N } )
        {
            auto result = mh.step( state, engine );
            OPENGEODE_EXCEPTION(
                result.decision == geode::MHDecision::Accepted
                    || result.decision == geode::MHDecision::Rejected,
                "[MH test] decision should be Accepted or Rejected." );

            // Log each step (optional: comment out if too verbose)
            //            geode::Logger::info( "Step: ", count,
            //                " move_type= ", static_cast< int >(
            //                result.move_type ), " decision= ", result.decision
            //                == geode::MHDecision::Accepted ? "Accepted"
            //                                                               :
            //                                                               "Rejected",
            //                " delta_log_energy = ", result.delta_log_energy,
            //                " log_accept = ", result.log_accept,
            //                " state_size = ", state.size() );

            // Keep track of accepted moves by type
            if( result.decision == geode::MHDecision::Accepted )
            {
                nb_accepted++;
                switch( result.move_type )
                {
                    case geode::Proposal< geode::Point2D >::Type::Birth:
                        accepted_birth++;
                        break;
                    case geode::Proposal< geode::Point2D >::Type::Death:
                        accepted_death++;
                        break;
                    case geode::Proposal< geode::Point2D >::Type::Change:
                        accepted_change++;
                        break;
                    default:
                        break;
                }
            }
            // should be change... only pone group here
            stat_sum += state.nb_objects_in_group( group_id );

            if( count % 1000 == 0 )
            {
                geode::Logger::info( "Progress: ", count, "  ", N,
                    " Mean objects =  ",
                    static_cast< double >( stat_sum )
                        / static_cast< double >( count ),
                    " nb accepted = ", nb_accepted, " Accepted(B/D/C) = ",
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
    }

} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();

        geode::Point2D min_point{ { 0., 0. } };
        geode::Point2D max_point{ { 10., 10. } };

        geode::BoundingBox2D box;
        box.add_point( min_point );
        box.add_point( max_point );

        geode::uuid group_id;
        geode::UniformPointConfigurationSampler< 2 > sampler( box, group_id );
        double birth_prob = 0.3;
        double death_prob = 0.1;
        auto kernel = geode::create_birth_death_change_kernel< geode::Point2D >(
            sampler, birth_prob, death_prob );

        geode::GibbsEnergy< geode::Point2D > poisson_energy;

        // Add intensity term
        poisson_energy.add_energy_term(
            std::make_unique< geode::IntensityTerm< geode::Point2D > >(
                0.5, group_id ) );

        geode::MetropolisHastings< geode::Point2D > mh(
            poisson_energy, std::move( kernel ) );

        test_steps( mh, group_id );
        test_beta_setter( mh );
        test_acceptance_prob_helper();

        geode::Logger::info( "MH TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}