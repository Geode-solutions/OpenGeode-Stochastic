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
#include <geode/stochastic/configuration/configuration.hpp>
#include <geode/stochastic/sampling/direct/marked_object_sampler/uniform_marked_point_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/moves.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/proposal_kernel.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

namespace
{
    geode::Configuration< geode::Point2D > create_configuration()
    {
        geode::Point2D p1{ { 0., 0. } };
        geode::MarkedObject< geode::Point2D > mp1{ std::move( p1 ) };
        geode::Point2D p2{ { 1., 1. } };
        geode::MarkedObject< geode::Point2D > mp2{ std::move( p2 ) };

        geode::Configuration< geode::Point2D > pattern;
        pattern.add_object( std::move( mp1 ) );
        pattern.add_object( std::move( mp2 ) );

        return pattern;
    }
} // namespace

void test_proposal_kernel()
{
    geode::Configuration< geode::Point2D > config = create_configuration();
    geode::Point2D min_point{ { 0., 0. } };
    geode::Point2D max_point{ { 10., 100. } };

    geode::BoundingBox2D box;
    box.add_point( min_point );
    box.add_point( max_point );

    geode::UniformMarkedPointSampler< 2 > sampler( box, std::nullopt );

    // Create classical birth-death-change kernel
    auto kernel = geode::create_birth_death_change_kernel< geode::Point2D >(
        sampler, 0.5 );

    geode::RandomEngine engine;

    for( int i = 0; i < 50; ++i )
    {
        auto prop = kernel->propose( config, engine );

        // Validate proposal type
        switch( prop.type )
        {
            case geode::Proposal< geode::Point2D >::Type::Birth:
                OPENGEODE_EXCEPTION( prop.new_object.has_value(),
                    "[test proposal] Birth must provide new_object." );
                OPENGEODE_EXCEPTION( !prop.index.has_value(),
                    "[test proposal] Birth should not provide index." );
                break;

            case geode::Proposal< geode::Point2D >::Type::Death:
                OPENGEODE_EXCEPTION( !prop.new_object.has_value(),
                    "[test proposal] Death should not provide new_object." );
                OPENGEODE_EXCEPTION( prop.index.has_value(),
                    "[test proposal] Death must provide index." );
                OPENGEODE_EXCEPTION( prop.index.value() < config.size(),
                    "[test proposal] Death index out of bounds." );
                break;

            case geode::Proposal< geode::Point2D >::Type::Change:
                OPENGEODE_EXCEPTION( prop.new_object.has_value(),
                    "[test proposal] Change must provide new_object." );
                OPENGEODE_EXCEPTION( prop.index.has_value(),
                    "[test proposal] Change must provide index." );
                OPENGEODE_EXCEPTION( prop.index.value() < config.size(),
                    "[test proposal] Change index out of bounds." );
                break;
        }

        // Log probabilities must be finite
        OPENGEODE_EXCEPTION( std::isfinite( prop.log_forward_prob ),
            "[test proposal] Forward probability is not finite." );
        OPENGEODE_EXCEPTION( std::isfinite( prop.log_backward_prob ),
            "[test proposal] Backward probability is not finite." );
    }
}

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