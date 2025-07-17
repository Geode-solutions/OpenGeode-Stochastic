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

#include <geode/stochastic/common.hpp>

#include <optional>

namespace geode
{
    template < typename Type >
    struct UniformClosed
    {
        UniformClosed() = default;
        bool is_valid() const;

        Type min_value{ static_cast< Type >( 0 ) };
        Type max_value{ static_cast< Type >( 1 ) };
    };

    template < typename Type >
    struct UniformClosedOpen
    {
        UniformClosedOpen() = default;
        bool is_valid() const;

        Type min_value{ static_cast< Type >( 0 ) };
        Type max_value{ static_cast< Type >( 1 ) };
    };

    struct opengeode_stochastic_stochastic_api Gaussian
    {
        Gaussian() = default;
        bool is_valid() const;

        double mean{ 0. };
        double standard_deviation{ 1. };
    };

    struct opengeode_stochastic_stochastic_api TruncatedGaussian
    {
        TruncatedGaussian() = default;

        bool is_valid() const;

        double mean{ 0. };
        double standard_deviation{ 1. };

        std::optional< double > min_value;
        std::optional< double > max_value;
    };

} // namespace geode