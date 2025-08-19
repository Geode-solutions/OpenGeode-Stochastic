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

#include <geode/stochastic/sampling/random_engine.hpp>

namespace geode
{
    template < typename Geometry >
    class MarkedObjectSampler
    {
    public:
        using MarkedObj = MarkedObject< Geometry >;

        virtual ~MarkedObjectSampler() = default;

        // Sample a new MarkedObject given an optional RNG
        virtual MarkedObj sample( RandomEngine& engine ) const = 0;
        virtual MarkedObject< Geometry > change(
            const MarkedObject< Geometry >& obj,
            RandomEngine& engine ) const = 0;

        // Log-probability density of a given object under this sampler
        virtual double log_pdf( const MarkedObj& obj ) const = 0;
    };

    //    template <typename Geometry>
    // class MultiTypeSampler
    //{
    // public:
    //    void add_sampler(std::unique_ptr<MarkedObjectSampler<Geometry>>
    //    sampler)
    //    {
    //        samplers_.push_back(std::move(sampler));
    //    }
    //
    //    const MarkedObjectSampler<Geometry>& get_sampler(std::size_t type)
    //    const
    //    {
    //        return *samplers_.at(type);
    //    }
    //
    //    std::size_t num_types() const { return samplers_.size(); }
    //
    // private:
    //    std::vector<std::unique_ptr<MarkedObjectSampler<Geometry>>> samplers_;
    //};
} // namespace geode