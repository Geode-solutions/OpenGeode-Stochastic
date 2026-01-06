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

#pragma once
#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace geode
{
    enum class MoveType
    {
        Invalid,
        Birth,
        Death,
        Change
    };

    struct ProposalProbabilities
    {
        double log_forward_prob{ 0. };
        double log_backward_prob{ 0. };
        std::string string() const
        {
            return absl::StrCat(
                "\t Proposal Probabilities: (log) forward probability: ",
                log_forward_prob,
                " (log) backward probability: ", log_backward_prob );
        }

        double transition_probability() const
        {
            OPENGEODE_EXCEPTION( std::isfinite( log_forward_prob )
                                     && std::isfinite( log_backward_prob ),
                "[Proposal Probabilities] Non-finite proposal "
                "log-probabilities" );
            return log_backward_prob - log_forward_prob;
        }
    };

    template < typename ObjectType >
    struct MoveResult
    {
        MoveType type{ MoveType::Invalid };
        std::optional< ObjectType > new_object;
        std::optional< index_t > old_object_id;

        ProposalProbabilities proposal_probabilities;

        std::string type_string() const
        {
            switch( type )
            {
                case MoveType::Invalid:
                    return "Invalid";
                case MoveType::Birth:
                    return "Birth";
                case MoveType::Death:
                    return "Death";
                case MoveType::Change:
                    return "Change";
                default:
                    throw OpenGeodeException(
                        "[MoveResult] -Move result type should always be "
                        "defined." );
                    return "Unknown";
            }
        }
        std::string string() const
        {
            auto message = absl::StrCat( "Move result type: ", type_string() );
            if( new_object )
            {
                absl::StrAppend( &message, "\t - New proposed object: ",
                    geode::object_bounding_box( new_object->first ).string() );
            }
            if( old_object_id )
            {
                absl::StrAppend(
                    &message, "\t - Old object id: ", old_object_id );
            }
            absl::StrAppend(
                &message, "\t - ", proposal_probabilities.string() );
            return message;
        }
    };

    template < typename ObjectType >
    class Move
    {
    public:
        Move( const ObjectSetSampler< ObjectType > &sampler,
            double proportion_weight )
            : sampler_( sampler ), proportion_weight_{ proportion_weight }
        {
            OPENGEODE_EXCEPTION( proportion_weight_ > 0.,
                "[Move] - the weight factor for a move should be in higher "
                "than 0. (here = ",
                proportion_weight_, ")" );
            initialize_probability( 1. );
        }
        virtual ~Move() = default;

        virtual MoveResult< ObjectType > propose_move(
            const ObjectSet< ObjectType > &set,
            RandomEngine &engine ) const = 0;

        double proportion_weight() const
        {
            return proportion_weight_;
        }

        virtual void initialize_probability( double probability )
        {
            geode_unused( probability );
        }

        virtual std::string string() const = 0;

    protected:
        std::optional< geode::index_t > draw_a_free_sample_id(
            const ObjectSet< ObjectType > &set, RandomEngine &engine ) const
        {
            const auto max_obj_id = set.nb_free_objects();
            if( max_obj_id == 0 )
            {
                return std::nullopt;
            }
            geode::UniformClosed< index_t > uniform_closed_index_t;
            uniform_closed_index_t.min_value = 0;
            uniform_closed_index_t.max_value = max_obj_id - 1;
            return engine.sample_uniform( uniform_closed_index_t );
        }

    protected:
        const ObjectSetSampler< ObjectType > &sampler_;
        double proportion_weight_{ 1.0 };
    };

    template < typename ObjectType >
    class BirthDeathMove : public Move< ObjectType >
    {
    public:
        BirthDeathMove( const ObjectSetSampler< ObjectType > &sampler,
            double probability,
            double birth_ratio )
            : Move< ObjectType >( sampler, probability ),
              birth_ratio_( birth_ratio )
        {
        }

        MoveResult< ObjectType > propose_move(
            const ObjectSet< ObjectType > &set,
            RandomEngine &engine ) const override
        {
            if( engine.sample_bernoulli( birth_ratio_ ) )
            {
                return propose_birth_move( set, engine );
            }
            return propose_death_move( set, engine );
        }

        void initialize_probability( double probability ) override
        {
            OPENGEODE_EXCEPTION( birth_ratio_ > 0. && birth_ratio_ < 1.,
                "[BirthDeathMove]-the ratio of birth over mover should be in "
                "]0,1[. (here = ",
                birth_ratio_, ")" );
            log_p_birth_ = std::log( probability * birth_ratio_ );
            log_p_death_ = std::log( probability * ( 1.0 - birth_ratio_ ) );
        }
        std::string string() const override
        {
            return absl::StrCat( "Birth and Death Move (proportion weight: ",
                this->proportion_weight_, " -- birth ratio:", birth_ratio_,
                " log_p_birth: ", log_p_birth_, " log_p_death_: ", log_p_death_,
                ")" );
        }

    private:
        MoveResult< ObjectType > propose_birth_move(
            const ObjectSet< ObjectType > &set, RandomEngine &engine ) const
        {
            MoveResult< ObjectType > birth;
            birth.type = MoveType::Birth;
            birth.new_object = this->sampler_.sample( engine );
            if( !birth.new_object.has_value() )
            {
                return birth;
            }
            auto &new_obj = birth.new_object.value();
            birth.proposal_probabilities.log_forward_prob =
                log_p_birth_ + this->sampler_.log_pdf( new_obj );
            birth.proposal_probabilities.log_backward_prob =
                log_p_death_ - std::log( set.nb_objects() + 1.0 );
            return birth;
        }

        MoveResult< ObjectType > propose_death_move(
            const ObjectSet< ObjectType > &set, RandomEngine &engine ) const
        {
            MoveResult< ObjectType > death;
            death.old_object_id = this->draw_a_free_sample_id( set, engine );
            if( !death.old_object_id.has_value() )
            {
                return death;
            }
            const auto cur_object_id = death.old_object_id.value();
            death.type = MoveType::Death;
            death.proposal_probabilities.log_forward_prob =
                log_p_death_ - std::log( set.nb_free_objects() );
            death.proposal_probabilities.log_backward_prob =
                log_p_birth_
                + this->sampler_.log_pdf(
                    set.get_free_object( cur_object_id ) );
            return death;
        }

    private:
        double birth_ratio_{ 0.5 };
        double log_p_birth_{ 0. };
        double log_p_death_{ 0. };
    };

    template < typename ObjectType >
    class ChangeMove : public Move< ObjectType >
    {
    public:
        ChangeMove(
            const ObjectSetSampler< ObjectType > &sampler, double probability )
            : Move< ObjectType >( sampler, probability )
        {
        }

        MoveResult< ObjectType > propose_move(
            const ObjectSet< ObjectType > &set,
            RandomEngine &engine ) const override
        {
            MoveResult< ObjectType > change;
            change.old_object_id = this->draw_a_free_sample_id( set, engine );
            if( !change.old_object_id.has_value() )
            {
                return change;
            }
            change.type = MoveType::Change;
            const auto &object_to_change =
                set.get_free_object( change.old_object_id.value() );
            change.new_object =
                this->sampler_.change( object_to_change, engine );
            change.proposal_probabilities.log_forward_prob =
                this->sampler_.log_pdf( change.new_object.value() );
            change.proposal_probabilities.log_backward_prob =
                this->sampler_.log_pdf( object_to_change );
            return change;
        }
        std::string string() const override
        {
            return absl::StrCat( "Change Move (proportion weight: ",
                this->proportion_weight_, ")" );
        }
    };
} // namespace geode