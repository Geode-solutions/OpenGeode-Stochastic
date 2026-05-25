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

#pragma once
#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <memory>
#include <vector>

namespace geode
{

    template < typename ObjectType >
    class StateDynamics
    {
    public:
        StateDynamics() = default;
        ~StateDynamics() = default;

        ObjectSets< ObjectType >& state()
        {
            return object_sets_;
        }
        const ObjectSets< ObjectType >& state() const
        {
            return object_sets_;
        }

        std::vector< std::unique_ptr< ObjectSetSampler< ObjectType > > >&
            samplers()
        {
            return samplers_;
        }
        const std::vector< std::unique_ptr< ObjectSetSampler< ObjectType > > >&
            samplers() const
        {
            return samplers_;
        }

        uuid add_set( const std::string& name )
        {
            return object_sets_.add_set( name );
        }

        void add_sampler(
            std::unique_ptr< ObjectSetSampler< ObjectType > > sampler )
        {
            samplers_.push_back( std::move( sampler ) );
        }

        ObjectSetSampler< ObjectType >& sampler( const index_t sampler_id )
        {
            OpenGeodeStochasticStochasticException::check_exception(
                sampler_id < samplers_.size(), nullptr,
                OpenGeodeException::TYPE::data,
                "[STATE DYNAMICS]: Sampler out of range." );
            return samplers_[sampler_id].get();
        }

    private:
        ObjectSets< ObjectType > object_sets_;
        std::vector< std::unique_ptr< ObjectSetSampler< ObjectType > > >
            samplers_;
    };

} // namespace geode