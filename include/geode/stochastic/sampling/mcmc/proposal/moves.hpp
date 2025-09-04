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
#include <geode/stochastic/configuration/configuration.hpp>
#include <geode/stochastic/sampling/direct/configuration_sampler/configuration_sampler.hpp>

namespace geode
{
    template < typename Object >
    struct Proposal
    {
        enum class Type
        {
            Invalid,
            Birth,
            Death,
            Change
        };

        Type type{ Type::Invalid };
        std::optional< std::pair< Object, GroupId > >
            new_object; // for birth/change
        std::optional< ObjectId > old_object_id; // for death/change
        double log_forward_prob{ 0. };
        double log_backward_prob{ 0. };
    };

    // Move does not hold the sampler... should it?
    template < typename Object >
    class Move
    {
    public:
        Move(
            const ConfigurationSampler< Object >& sampler, double probability )
            : sampler_( sampler ), p_move_{ probability }
        {
        }
        virtual ~Move() = default;

        virtual Proposal< Object > propose_move(
            const Configuration< Object >& current,
            RandomEngine& engine ) const = 0;

        double probability() const
        {
            return p_move_;
        }

    protected:
        const ConfigurationSampler< Object >& sampler_;
        double p_move_{ 1.0 };
    };

    template < typename Object >
    class BirthDeathMove : public Move< Object >
    {
    public:
        BirthDeathMove( const ConfigurationSampler< Object >& sampler,
            double probability,
            double birth_ratio )
            : Move< Object >( sampler, probability ),
              birth_ratio_( birth_ratio )
        {
            OPENGEODE_EXCEPTION( birth_ratio_ > 0. && birth_ratio_ < 1.,
                "[BirthDeathMove]-the ratio of birth over mover should be in "
                "]0,1[. (here = ",
                birth_ratio_, ")" );
            log_p_birth_ = std::log( this->p_move_ * birth_ratio );
            log_p_death_ = std::log( this->p_move_ * ( 1.0 - birth_ratio ) );
        }

        Proposal< Object > propose_move( const Configuration< Object >& current,
            RandomEngine& engine ) const override
        {
            if( engine.sample_bernoulli( birth_ratio_ ) )
            {
                return propose_birth_move( current, engine );
            }
            return propose_death_move( current, engine );
        }

    private:
        Proposal< Object > propose_birth_move(
            const Configuration< Object >& current, RandomEngine& engine ) const
        {
            Proposal< Object > birth;
            birth.type = Proposal< Object >::Type::Birth;
            birth.new_object = this->sampler_.sample( engine );
            if( !birth.new_object.has_value() )
            {
                return birth;
            }
            auto& [new_obj, group_id] = birth.new_object.value();
            birth.log_forward_prob =
                log_p_birth_ + this->sampler_.log_pdf( new_obj );
            birth.log_backward_prob =
                log_p_death_
                - std::log( current.nb_objects_in_group( group_id ) + 1.0 );
            return birth;
        }

        Proposal< Object > propose_death_move(
            const Configuration< Object >& current, RandomEngine& engine ) const
        {
            Proposal< Object > death;
            death.old_object_id = this->sampler_.sample_id( current, engine );
            if( !death.old_object_id.has_value() )
            {
                return death;
            }
            const auto& cur_object_id = death.old_object_id.value();
            death.type = Proposal< Object >::Type::Death;
            death.log_forward_prob = log_p_death_
                                     - std::log( current.nb_objects_in_group(
                                         cur_object_id.group ) );
            death.log_backward_prob =
                log_p_birth_
                + this->sampler_.log_pdf( current.get_object( cur_object_id ) );
            return death;
        }

    private:
        double birth_ratio_{ 0.5 };
        double log_p_birth_{ 0. };
        double log_p_death_{ 0. };
    };

    template < typename Object >
    class ChangeMove : public Move< Object >
    {
    public:
        ChangeMove(
            const ConfigurationSampler< Object >& sampler, double probability )
            : Move< Object >( sampler, probability )
        {
        }

        Proposal< Object > propose_move( const Configuration< Object >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Object > change;
            change.old_object_id = this->sampler_.sample_id( current, engine );
            if( !change.old_object_id.has_value() )
            {
                return change;
            }
            change.type = Proposal< Object >::Type::Change;
            const auto& object_to_change =
                current.get_object( change.old_object_id.value() );
            change.new_object =
                this->sampler_.change( object_to_change, engine );
            change.log_forward_prob =
                this->sampler_.log_pdf( change.new_object.value().first );
            change.log_backward_prob =
                this->sampler_.log_pdf( object_to_change );
            return change;
        }
    };
} // namespace geode