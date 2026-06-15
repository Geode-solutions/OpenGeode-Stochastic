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

#include <geode/basic/uuid.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/distance.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace geode
{
    template < typename Type >
    class PairwiseInteraction
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( PairwiseInteraction );

    public:
        explicit PairwiseInteraction() = default;
        virtual ~PairwiseInteraction() = default;

        [[nodiscard]] double evaluate( const ObjectRef< Type >& object_a,
            const ObjectRef< Type >& object_b ) const
        {
            return compute( object_a, object_b );
        }

        [[nodiscard]] virtual double
            neighborhood_searching_distance() const = 0;

    protected:
        [[nodiscard]] virtual double compute( const ObjectRef< Type >& object_a,
            const ObjectRef< Type >& object_b ) const = 0;
    };

} // namespace geode