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

#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>

namespace
{
    void register_interaction(
        const std::pair< geode::uuid, geode::uuid >& interaction,
        std::vector< geode::uuid >& interacting_set_ids,
        absl::flat_hash_map< geode::uuid, std::vector< geode::uuid > >&
            adjacency_map )
    {
        const auto& [set1, set2] = interaction;

        adjacency_map[set1].push_back( set2 );
        adjacency_map[set2].push_back( set1 );

        interacting_set_ids.push_back( set1 );
        interacting_set_ids.push_back( set2 );
    }

    void sort_unique_shrink( std::vector< geode::uuid >& values )
    {
        absl::c_sort( values );
        values.erase(
            std::unique( values.begin(), values.end() ), values.end() );
        values.shrink_to_fit();
    }

    void deduplicate_adjacency_map(
        absl::flat_hash_map< geode::uuid, std::vector< geode::uuid > >&
            adjacency_map )
    {
        for( auto& [set_uuid, adjacent_set_uuids] : adjacency_map )
        {
            geode_unused( set_uuid );
            sort_unique_shrink( adjacent_set_uuids );
        }
    }
} // namespace

namespace geode
{
    std::pair< std::vector< geode::uuid >,
        absl::flat_hash_map< uuid, std::vector< uuid > > >
        pairwise_builder_initialize_interactions_helper(
            const std::vector< std::pair< uuid, uuid > >& interacting_sets )
    {
        std::vector< geode::uuid > interacting_set_ids;
        absl::flat_hash_map< uuid, std::vector< uuid > >
            objectset_adjacency_map;

        interacting_set_ids.reserve( 2 * interacting_sets.size() );
        objectset_adjacency_map.reserve( 2 * interacting_sets.size() );

        for( const auto& interaction : interacting_sets )
        {
            register_interaction(
                interaction, interacting_set_ids, objectset_adjacency_map );
        }

        sort_unique_shrink( interacting_set_ids );
        deduplicate_adjacency_map( objectset_adjacency_map );

        return { std::move( interacting_set_ids ),
            std::move( objectset_adjacency_map ) };
    }
} // namespace geode
