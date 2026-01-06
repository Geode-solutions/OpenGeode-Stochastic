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
    public:
        enum struct SCOPE
        {
            all_set,
            same_set, // only within the same subset
            different_set // only across subsets
        };

        explicit PairwiseInteraction() = default;
        explicit PairwiseInteraction( SCOPE scope ) : scope_( scope ) {}
        virtual ~PairwiseInteraction() = default;

        double evaluate( const ObjectRef< Type >& object_a,
            const ObjectRef< Type >& object_b ) const
        {
            if( scope_ == SCOPE::same_set
                && object_a.set_id != object_b.set_id )
            {
                return 0.0;
            }
            if( scope_ == SCOPE::different_set
                && object_a.set_id == object_b.set_id )
            {
                return 0.0;
            }
            return compute( object_a, object_b );
        }

        virtual double neighborhood_searching_distance() const = 0;

    protected:
        virtual double compute( const ObjectRef< Type >& object_a,
            const ObjectRef< Type >& object_b ) const = 0;

    private:
        SCOPE scope_{ SCOPE::all_set };
    };

    /*!
     * EuclideanCutoffInteraction
     * A pairwise interaction that returns 1 if the Euclidean distance
     * between objects is within a cutoff radius, otherwise 0.
     */
    template < typename Type >
    class EuclideanCutoffInteraction : public PairwiseInteraction< Type >
    {
    public:
        explicit EuclideanCutoffInteraction( double cutoff_distance );
        EuclideanCutoffInteraction( double cutoff_distance,
            typename PairwiseInteraction< Type >::SCOPE scope );

        double neighborhood_searching_distance() const override;

    protected:
        double compute( const ObjectRef< Type >& object_a,
            const ObjectRef< Type >& object_b ) const override;

    private:
        double cutoff_distance_{ GLOBAL_EPSILON };
    };
} // namespace geode