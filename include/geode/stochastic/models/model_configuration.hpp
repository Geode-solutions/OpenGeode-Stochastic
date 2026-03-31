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

#include <absl/container/btree_map.h>
#include <absl/strings/str_join.h>

#include <geode/basic/uuid.hpp>

#include <string>
#include <variant>
#include <vector>

#include <geode/stochastic/models/energy_term_collection.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>
namespace geode
{

    struct ModelConfig
    {
        std::vector< EnergyTermConfig > terms;
    };

    template < typename ObjectType >
    EnergyTermCollection< ObjectType > build_energy_term_collection(
        const ModelConfig& config,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        EnergyTermCollection< ObjectType > collection;

        for( const auto& term_cfg : config.terms )
        {
            auto term_id = collection.add_energy_term(
                build_term< ObjectType >( term_cfg, object_sets, domain ) );
        }

        return collection;
    }
} // namespace geode