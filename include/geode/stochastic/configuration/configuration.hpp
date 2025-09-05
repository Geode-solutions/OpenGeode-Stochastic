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

#include <vector>

#include <geode/basic/uuid.hpp>
#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

// #include <geode/stochastic/configuration/neighbors_object_ids.hpp>
namespace geode
{
    struct GroupDescriptor
    {
        geode::uuid unique_id;
        bool operator==( GroupDescriptor const& other ) const noexcept
        {
            {
                return unique_id == other.unique_id;
            }
        }
    };

    struct ObjectId
    {
        index_t object;
        uuid group;
        bool operator==( const ObjectId& other ) const noexcept
        {
            return object == other.object && group == other.group;
        }
    };
} // namespace geode
// Hash support for unordered_map
namespace std
{
    // template <>
    // struct hash< geode::uuid >
    //{
    //     std::size_t operator()( geode::uuid const& g ) const noexcept
    //     {
    //         return std::hash< geode::index_t >()( g.unique_id );
    //     }
    // };

    // template <>
    // struct hash< geode::ObjectId >
    //{
    //     std::size_t operator()( geode::ObjectId const& id ) const noexcept
    //     {
    //         return ( std::hash< geode::index_t >()( id.group.value ) << 1 )
    //                ^ std::hash< geode::index_t >()( id.object );
    //     }
    // };
} // namespace std
namespace geode
{
    template < typename Object >
    class Configuration
    {
    public:
        const std::vector< Object >& get_group( const uuid& group_id ) const;
        const Object& get_object( const ObjectId& object_id ) const;
        std::vector< ObjectId > get_all_object() const;

        index_t nb_groups() const;
        index_t nb_objects_in_group( const uuid& group_id ) const;
        index_t nb_objects() const;

        void add_group( const uuid& group_id );
        ObjectId add_object( Object&& object, const uuid& group_id );
        void update_object( const ObjectId& object_id, Object&& object );
        void remove_object( const ObjectId& object_id );

        // Object neighbor search by ObjectId (always excludes self)
        std::vector< ObjectId > neighbors(
            const ObjectId& object_id, double searching_distance ) const;
        // Object neighbor search by arbitrary object (return self if in the
        // configuration)
        std::vector< ObjectId > neighbors(
            const Object& object, double searching_distance ) const;

    private:
        std::vector< Object >& get_group( const uuid& group_id );

    private:
        std::unordered_map< uuid, std::vector< Object > > groups_;
        // ObjectIndexRTree< 2 > tree_;
    };
} // namespace geode
