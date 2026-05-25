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

#include <string>
#include <variant>
#include <vector>

#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions_config.hpp>
#include <geode/stochastic/spatial/single_object_features/single_object_feature_config.hpp>

namespace geode
{
    struct SingleObjectTermConfig
    {
        std::string term_name;
        std::vector< std::string > object_set_names;
        double lambda;

        SingleObjectFeatureConfig object_feature;
    };

    struct PairwiseTermConfig
    {
        std::string term_name;
        std::vector< std::pair< std::string, std::string > >
            object_set_names_interactions;
        double gamma;

        PairwiseInteractionConfig interaction_config;
    };

    using EnergyTermConfig = std::
        variant< std::monostate, SingleObjectTermConfig, PairwiseTermConfig >;

} // namespace geode