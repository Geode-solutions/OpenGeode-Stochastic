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
#include <optional>

namespace geode
{
    template < typename Type >
    struct IntervalLimit
    {
        IntervalLimit() = default;

        Type value{ 0 };
        bool is_included{ true };
    };

    template < typename Type >
    struct Uniform
    {
        Uniform() = default;
        Uniform( Type min_val, Type max_val )
            : min{ IntervalLimit< Type >{ min_val, true } },
              max{ IntervalLimit< Type >{ max_val, true } }
        {
        }
        IntervalLimit< Type > min;
        IntervalLimit< Type > max;
    };

    struct Gaussian
    {
        Gaussian() = default;
        Gaussian( double mean_val, double standard_deviation_val )
            : mean{ mean_val }, standard_deviation{ standard_deviation_val }
        {
        }
        Gaussian( double mean_val,
            double standard_deviation_val,
            double min_val,
            double max_val )
            : mean{ mean_val },
              standard_deviation{ standard_deviation_val },
              min{ min_val },
              max{ max_val }
        {
        }

        double mean{ 0. };
        double standard_deviation{ 1. };

        std::optional< double > min;
        std::optional< double > max;
    };

} // namespace geode