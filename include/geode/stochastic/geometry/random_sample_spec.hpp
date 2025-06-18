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

#include <geode/stochastic/geometry/common.hpp>

#include <geode/basic/pimpl.hpp>

#include <variant>

namespace geode
{
    /*!
     * Spec to draw a int value in a Uniform Distribution. Describe a
     * closed interval [min,max]
     */
    struct UniformInt
    {
        int min;
        int max;
    };

    /*!
     * Spec to draw a double value in a Uniform Distribution. Describe a
     * semi closed interval [min,max[
     */
    struct UniformDouble
    {
        double min;
        double max;
    };

    /*!
     * Spec to draw a double value in a Gaussian Distribution defined by its
     * mean and standard deviation.
     */
    struct Gaussian
    {
        double mean;
        double stddev;
    };

    /*!
     * Spec to draw a sucess/fail of an experiement. the experiment is defined
     * by its probability of success.
     */
    struct Bernoulli
    {
        double probability;
    };

    //    struct PointInDisk
    //    {
    //        double radius;
    //    };
    //
    //    struct PointInBox2D
    //    {
    //        geode::BoundingBox2D box;
    //    };

    using RandomDoubleSpec = std::variant< UniformDouble, Gaussian >;
    //        PointInDisk,
    //        PointInBox2D >;
} // namespace geode