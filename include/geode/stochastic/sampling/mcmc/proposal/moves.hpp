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
#include <geode/stochastic/sampling/direct/marked_object_sampler/marked_object_sampler.hpp>
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

        Type type{ Undefined };
        std::optional< MarkedObject< Geometry > >
            new_object; // for birth/change
        std::optional< index_t > index; // for death/change
        double log_forward_prob{ 1. };
        double log_backward_prob{ 1. };
    };

    template < typename Geometry >
    class Move
    {
    public:
        Move( double probability ) : probability_{ probability } {}
        virtual ~Move() = default;

        // Propose a new move given the current configuration and RNG
        virtual Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const = 0;

        double probability() const
        {
            return probability_;
        }

    protected:
        double probability_{ 1.0 };
    };

    template < typename Geometry >
    class BirthMove : public Move< Geometry >
    {
    public:
        BirthMove(
            double probability, const MarkedObjectSampler< Geometry >& sampler )
            : Move< Geometry >( probability ), sampler_( sampler )
        {
        }

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > p;
            p.type = Proposal< Geometry >::Type::Birth;
            p.new_object = sampler_.sample( engine );
            p.index = std::nullopt;
            p.log_forward_prob = sampler_.log_pdf( *p.new_object );
            p.log_backward_prob = -std::log( current.size() + 1.0 );
            return p;
        }

    private:
        const MarkedObjectSampler< Geometry >& sampler_;
    };

    template < typename Geometry >
    class DeathMove : public Move< Geometry >
    {
    public:
        DeathMove( double probability ) : Move< Geometry >( probability ) {}

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > p;
            if( current.size() == 0 )
                return p;

            geode::UniformClosed< index_t > uniform_closed_index_t;
            uniform_closed_index_t.min_value = 0;
            uniform_closed_index_t.max_value = current.size() - 1;
            index_t idx = engine.sample_uniform( uniform_closed_index_t );

            p.type = Proposal< Geometry >::Type::Death;
            p.index = idx;
            p.log_forward_prob = -std::log( current.size() );
            p.log_backward_prob = 0; // Could depend on birth sampler if needed
            p.new_object = std::nullopt;
            return p;
        }
    };

    template < typename Geometry >
    class ChangeMove : public Move< Geometry >
    {
    public:
        ChangeMove(
            double probability, const MarkedObjectSampler< Geometry >& sampler )
            : Move< Geometry >( probability ), sampler_( sampler )
        {
        }

        Proposal< Geometry > propose_move(
            const Configuration< Geometry >& current,
            RandomEngine& engine ) const override
        {
            Proposal< Geometry > p;
            if( current.size() == 0 )
            {
                return p;
            }
            geode::UniformClosed< index_t > uniform_closed_index_t;
            uniform_closed_index_t.min_value = 0;
            uniform_closed_index_t.max_value = current.size() - 1;
            index_t idx = engine.sample_uniform( uniform_closed_index_t );

            p.type = Proposal< Geometry >::Type::Change;
            p.index = idx;
            p.new_object = sampler_.change( current.object( idx ), engine );
            p.log_forward_prob = sampler_.log_pdf( *p.new_object );
            p.log_backward_prob = sampler_.log_pdf( current.object( idx ) );
            return p;
        }

    private:
        const MarkedObjectSampler< Geometry >& sampler_;
    };
} // namespace geode