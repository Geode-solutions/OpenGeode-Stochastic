
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

#include <geode/stochastic/geometry/distributions.hpp>

#include <geode/stochastic/geometry/spatial_distributions.hpp>

#include <variant>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( Point );
    ALIAS_2D_AND_3D( Point );
    class RandomEngine;
} // namespace geode

namespace geode
{
    template < typename T >
    struct RandomSampler;

    using DoubleDistribution = std::variant< UniformClosed< double >,
        UniformClosedOpen< double >,
        Gaussian,
        TruncatedGaussian >;

    template <>
    struct RandomSampler< double >
    {
        static double sample(
            RandomEngine& engine, const DoubleDistribution& dist );
    };

    using SpatialDistribution2D =
        std::variant< UniformBox< 2 >, UniformBall< 2 > >;
    template <>
    struct RandomSampler< Point2D >
    {
        static Point2D sample(
            RandomEngine& engine, const SpatialDistribution2D& dist );
    };

    using SpatialDistribution3D =
        std::variant< UniformBox< 3 >, UniformBall< 3 > >;
    template <>
    struct RandomSampler< Point3D >
    {
        static Point3D sample(
            RandomEngine& engine, const SpatialDistribution3D& dist );
    };

} // namespace geode