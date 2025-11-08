
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

#include <geode/stochastic/sampling/distributions.hpp>

#include <variant>

namespace geode
{
    class RandomEngine;
} // namespace geode


namespace geode
{

    struct opengeode_stochastic_stochastic_api DoubleSampler
    {
        using Distribution = std::variant< UniformClosed< double >,
            UniformClosedOpen< double >,
            Gaussian,
            TruncatedGaussian >;

        struct DistributionDescription
        {
            std::string name{ "default_distribution" };
            DistributionType distribution_type;

            std::optional< double > min_value;
            std::optional< double > max_value;
            std::optional< double > mean;
            std::optional< double > standard_deviation;
        };

        static Distribution create_distribution(
            const DistributionDescription& desc );

        static double sample( RandomEngine& engine, const Distribution& dist );
    };

} // namespace geode