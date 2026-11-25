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

#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < typename ObjectType >
    class SingleObjectFeature
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( SingleObjectFeature );

    public:
        SingleObjectFeature() = default;
        virtual ~SingleObjectFeature() = default;

        /// Contribution of the object to the energy, wherever it is in the
        /// extended domain.
        [[nodiscard]] virtual double evaluate( const ObjectType& obj,
            const SpatialDomain< ObjectType::dim >& domain ) const = 0;

        /// Contribution of the object to the statistic, observed in the
        /// domain.
        [[nodiscard]] virtual double evaluate_inside_domain(
            const ObjectType& obj,
            const SpatialDomain< ObjectType::dim >& domain ) const = 0;
    };

    /// Counts the objects: 1 per object in the extended domain (energy), and
    /// the share of each object inside the domain (statistic), as defined by
    /// SpatialDomainChecker::count_inside_domain. For segments, this share is
    /// counted by endpoints and does not depend on the anchor.
    template < typename ObjectType >
    class ObjectInDomainFeature : public SingleObjectFeature< ObjectType >
    {
    public:
        [[nodiscard]] double evaluate( const ObjectType& /*obj*/,
            const SpatialDomain< ObjectType::dim >& /*domain*/ ) const override
        {
            return 1.;
        }

        [[nodiscard]] double evaluate_inside_domain( const ObjectType& obj,
            const SpatialDomain< ObjectType::dim >& domain ) const override
        {
            return SpatialDomainChecker< ObjectType >::count_inside_domain(
                domain, obj );
        }
    };

} // namespace geode