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

namespace geode
{
    template < index_t dimension >
    struct ObjectSamplerConfig< Point< dimension > >
    {
        // use to define the step for change move (move_ratio*domain volume)
        // NOLINTBEGING(*-magic-numbers)
        double move_ratio{ 0.1 };
        // NOLINTEND(*-magic-numbers)
    };

    template < index_t dimension >
    class UniformPointSetSampler : public ObjectSetSampler< Point< dimension > >
    {
    public:
        UniformPointSetSampler( const SpatialDomain< dimension >& domain,
            const ObjectSamplerConfig< Point< dimension > >& config );

        [[nodiscard]] Point< dimension > sample(
            RandomEngine& engine ) const override;

        [[nodiscard]] Point< dimension > change( const Point< dimension >& obj,
            RandomEngine& engine ) const override;

    private:
        [[nodiscard]] double define_step_for_move( double ratio );
        [[nodiscard]] bool is_valid_object(
            const Point< dimension >& obj ) const override;

    private:
        const SpatialDomain< dimension >& domain_;
        double step_move_{ 0. };
    };

} // namespace geode