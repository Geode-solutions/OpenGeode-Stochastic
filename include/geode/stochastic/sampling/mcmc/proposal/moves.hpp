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
#include <geode/stochastic/sampling/mcmc/proposal/marked_object_sampler/marked_object_sampler.hpp>
namespace geode
{
    template < typename Geometry >
    struct Proposal
    {
        enum class Type
        {
            Undefined,
            Birth,
            Death,
            Change
        };

        Type type{ Type::Undefined };
        std::optional< MarkedObject< Geometry > >
            new_object; // for birth/change
        std::optional< index_t > index; // for death/change
        double log_forward_prob{ 0. };
        double log_backward_prob{ 0. };
    };

    template < typename Geometry >
    class Move
    {
    public:
        Move(
            const MarkedObjectSampler< Geometry >& sampler, double probability )
            : sampler_( sampler ), probability_{ probability }
        {
        }
        virtual ~Move() = default;

        virtual Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const = 0;

        double probability() const
        {
            return probability_;
        }

    protected:
        const MarkedObjectSampler< Geometry >& sampler_;
        double probability_{ 1.0 };
    };

    template < typename Geometry >
    class BirthMove : public Move< Geometry >
    {
    public:
        BirthMove(
            const MarkedObjectSampler< Geometry >& sampler, double probability )
            : Move< Geometry >( sampler, probability )
        {
        }

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > birth;
            birth.type = Proposal< Geometry >::Type::Birth;
            birth.new_object = this->sampler_.sample( engine );
            birth.log_forward_prob =
                this->sampler_.log_pdf( birth.new_object.value() );
            birth.log_backward_prob = -std::log( current.size() + 1.0 );
            return birth;
        }
    };

    template < typename Geometry >
    class DeathMove : public Move< Geometry >
    {
    public:
        DeathMove(
            const MarkedObjectSampler< Geometry >& sampler, double probability )
            : Move< Geometry >( sampler, probability )
        {
        }

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > death;
            death.index = this->sampler_.sample_id( current, engine );
            if( !death.index.has_value() )
            {
                return death;
            }
            death.type = Proposal< Geometry >::Type::Death;
            death.log_forward_prob = -std::log( current.size() );
            death.log_backward_prob =
                this->sampler_.log_pdf( current[death.index.value()] );
            return death;
        }
    };

    template < typename Geometry >
    class ChangeMove : public Move< Geometry >
    {
    public:
        ChangeMove(
            const MarkedObjectSampler< Geometry >& sampler, double probability )
            : Move< Geometry >( sampler, probability )
        {
        }

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > change;
            change.index = this->sampler_.sample_id( current, engine );
            if( !change.index.has_value() )
            {
                return change;
            }
            change.type = Proposal< Geometry >::Type::Change;
            const auto& object_to_change = current[change.index.value()];
            change.new_object =
                this->sampler_.change( object_to_change, engine );
            change.log_forward_prob =
                this->sampler_.log_pdf( change.new_object.value() );
            change.log_backward_prob =
                this->sampler_.log_pdf( object_to_change );
            return change;
        }
    };
} // namespace geode