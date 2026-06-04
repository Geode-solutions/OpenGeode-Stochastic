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

#include <absl/container/flat_hash_map.h>

#include <geode/basic/uuid.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/spatial/object_neighborhood.hpp>
#include <geode/stochastic/spatial/object_set.hpp>

namespace geode
{
    template < typename Type >
    struct ObjectRef
    {
        const Type& object;
        uuid set_id;
    };
} // namespace geode

namespace geode
{
    template < typename Type >
    class ObjectSets
    {
        OPENGEODE_DISABLE_COPY( ObjectSets );

    public:
        ObjectSets() noexcept = default;
        ~ObjectSets() noexcept = default;

        ObjectSets( ObjectSets&& ) noexcept = default;
        ObjectSets& operator=( ObjectSets&& ) noexcept = default;

        [[nodiscard]] const ObjectSet< Type >& get_set(
            const uuid& set_id ) const;
        [[nodiscard]] const Type& get_object( const ObjectId& object_id ) const;
        [[nodiscard]] std::vector< ObjectId > get_all_object() const;
        [[nodiscard]] std::vector< ObjectId > get_objects_in_set(
            const uuid& set_id ) const;

        [[nodiscard]] index_t nb_sets() const;
        [[nodiscard]] index_t nb_objects_in_set( const uuid& set_id ) const;
        [[nodiscard]] index_t nb_objects() const;

        uuid add_set( std::string_view name );
        ObjectId add_object( Type&& object, const uuid& set_id, bool fixed );
        void update_free_object( const ObjectId& object_id, Type&& object );
        void remove_free_object( const ObjectId& object_id );

        [[nodiscard]] std::vector< ObjectId > neighbors( const Type& object,
            const std::vector< uuid >& targeted_set_ids,
            double searching_distance,
            std::optional< ObjectId > exclude_id ) const;

        [[nodiscard]] std::string string() const;

        [[nodiscard]] uuid get_set_uuid( std::string_view set_name ) const;
        [[nodiscard]] std::vector< uuid > get_set_uuids(
            const std::vector< std::string >& set_names ) const;
        [[nodiscard]] std::vector< std::pair< uuid, uuid > > get_set_uuid_pairs(
            const std::vector< std::pair< std::string, std::string > >&
                set_name_pairs ) const;

    private:
        [[nodiscard]] ObjectSet< Type >& get_set( const uuid& set_id );

    private:
        absl::flat_hash_map< std::string, uuid > name_to_uuid_;
        absl::flat_hash_map< uuid, index_t > uuid_to_index_;
        std::vector< ObjectSet< Type > > sets_;

        ObjectNeighborhood< Type::dim > neighborhood_;
    };
} // namespace geode
