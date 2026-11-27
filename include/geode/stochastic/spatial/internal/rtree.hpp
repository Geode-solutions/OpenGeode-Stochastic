/*
 * This code has been copied from https://github.com/nushoin/RTree
 *
 * Original code by Antonin Guttman and Michael Stonebraker (algorithm),
 * Melinda Green (ANSI C port), Paul Brook, Greg Douglas (templated C++ port),
 * Yariv Barkan and Gero Mueller.
 *
 * The RTree library is released into the Public Domain where such a
 * disclaimer is recognized, and under the MIT License otherwise:
 *
 * Copyright (c) 2011-2024 Yariv Barkan
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
 */

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/common.hpp>

namespace geode
{
    namespace internal
    {
        /*!
         * Dynamic R-tree (Guttman, 1984) of data indexed by their bounding
         * box, with the quadratic split and the spherical volume of the
         * original implementation.
         * Removing an entry requires its bounding box at insertion time.
         * @tparam Data default constructible, copyable and equality comparable
         * identifier of the indexed objects.
         */
        template < typename Data, index_t dimension, index_t max_nodes = 8 >
        class RTree
        {
            static_assert( dimension == 2 || dimension == 3,
                "[RTree] Only 2D and 3D trees are supported" );
            static_assert(
                max_nodes >= 2, "[RTree] A node needs at least two entries" );

            static constexpr index_t MIN_NODES = max_nodes / 2;
            static constexpr index_t NOT_TAKEN = NO_ID;

            struct Rect
            {
                std::array< double, dimension > min;
                std::array< double, dimension > max;
            };

            struct Node;

            /// A leaf branch holds data, an internal branch holds a subtree.
            struct Branch
            {
                Rect rect;
                std::unique_ptr< Node > child;
                Data data{};
            };

            struct Node
            {
                [[nodiscard]] bool is_leaf() const
                {
                    return level == 0;
                }

                index_t count{ 0 };
                /// Leaves are at level 0, their parents at level 1, ...
                index_t level{ 0 };
                std::array< Branch, max_nodes > branches;
            };

            /// Branches of an overflowing node, and their split in two groups
            struct Partition
            {
                std::array< Branch, max_nodes + 1 > branches;
                std::array< index_t, max_nodes + 1 > group;
                std::array< index_t, 2 > count{ 0, 0 };
                std::array< Rect, 2 > cover;
                std::array< double, 2 > area{ 0., 0. };
            };

        public:
            RTree() : root_{ std::make_unique< Node >() } {}

            void insert( const BoundingBox< dimension >& box, const Data& data )
            {
                Branch branch;
                branch.rect = to_rect( box );
                branch.data = data;
                insert_branch( std::move( branch ), 0 );
            }

            /// Remove the entry of the given data inserted with the given box.
            /// @return false if no such entry was found
            bool remove( const BoundingBox< dimension >& box, const Data& data )
            {
                std::vector< std::unique_ptr< Node > > eliminated_nodes;
                if( !remove_rec(
                        to_rect( box ), data, *root_, eliminated_nodes ) )
                {
                    return false;
                }
                // Reinsert the entries of the underfull nodes, last eliminated
                // first
                for( auto node = eliminated_nodes.rbegin();
                    node != eliminated_nodes.rend(); ++node )
                {
                    for( const auto b : Range{ ( *node )->count } )
                    {
                        insert_branch( std::move( ( *node )->branches[b] ),
                            ( *node )->level );
                    }
                }
                // Remove a root left with a single subtree
                if( root_->count == 1 && !root_->is_leaf() )
                {
                    root_ = std::move( root_->branches[0].child );
                }
                return true;
            }

            /*!
             * Call the callback on the data of each entry whose box intersects
             * (or touches) the given box.
             * @param callback bool( const Data& ), returning false to stop the
             * search
             */
            template < typename Callback >
            void search(
                const BoundingBox< dimension >& box, Callback&& callback ) const
            {
                search_rec( *root_, to_rect( box ), callback );
            }

        private:
            static Rect to_rect( const BoundingBox< dimension >& box )
            {
                Rect rect;
                for( const auto d : LRange{ dimension } )
                {
                    rect.min[d] = box.min().value( d );
                    rect.max[d] = box.max().value( d );
                    assert( rect.min[d] <= rect.max[d] );
                }
                return rect;
            }

            static bool overlap( const Rect& rect_a, const Rect& rect_b )
            {
                for( const auto d : LRange{ dimension } )
                {
                    if( rect_a.min[d] > rect_b.max[d]
                        || rect_b.min[d] > rect_a.max[d] )
                    {
                        return false;
                    }
                }
                return true;
            }

            static Rect combine( const Rect& rect_a, const Rect& rect_b )
            {
                Rect rect;
                for( const auto d : LRange{ dimension } )
                {
                    rect.min[d] = std::min( rect_a.min[d], rect_b.min[d] );
                    rect.max[d] = std::max( rect_a.max[d], rect_b.max[d] );
                }
                return rect;
            }

            /// Volume of the sphere circumscribed to the rectangle, which
            /// discriminates the splits better than the rectangle volume
            static double volume( const Rect& rect )
            {
                double sum_of_squares{ 0. };
                for( const auto d : LRange{ dimension } )
                {
                    const auto half_extent =
                        ( rect.max[d] - rect.min[d] ) * 0.5;
                    sum_of_squares += half_extent * half_extent;
                }
                const auto radius = std::sqrt( sum_of_squares );
                if constexpr( dimension == 2 )
                {
                    return radius * radius * UNIT_SPHERE_VOLUME;
                }
                else
                {
                    return radius * radius * radius * UNIT_SPHERE_VOLUME;
                }
            }

            static Rect node_cover( const Node& node )
            {
                auto rect = node.branches[0].rect;
                for( const auto b : Range{ 1, node.count } )
                {
                    rect = combine( rect, node.branches[b].rect );
                }
                return rect;
            }

            /// Insert a branch at the given level: 0 for data, the level of
            /// its subtree root otherwise. Grow the tree if the root splits.
            void insert_branch( Branch branch, index_t level )
            {
                assert( level <= root_->level );
                auto new_node =
                    insert_rec( std::move( branch ), *root_, level );
                if( !new_node )
                {
                    return;
                }
                auto new_root = std::make_unique< Node >();
                new_root->level = root_->level + 1;
                Branch old_root_branch;
                old_root_branch.rect = node_cover( *root_ );
                old_root_branch.child = std::move( root_ );
                add_branch( std::move( old_root_branch ), *new_root );
                Branch new_node_branch;
                new_node_branch.rect = node_cover( *new_node );
                new_node_branch.child = std::move( new_node );
                add_branch( std::move( new_node_branch ), *new_root );
                root_ = std::move( new_root );
            }

            /// @return the node created if the given node was split
            std::unique_ptr< Node > insert_rec(
                Branch branch, Node& node, index_t level )
            {
                assert( level <= node.level );
                if( node.level == level )
                {
                    return add_branch( std::move( branch ), node );
                }
                const auto index = pick_branch( branch.rect, node );
                auto& picked = node.branches[index];
                const auto rect = branch.rect;
                auto split_node =
                    insert_rec( std::move( branch ), *picked.child, level );
                if( !split_node )
                {
                    picked.rect = combine( rect, picked.rect );
                    return nullptr;
                }
                // The picked child was split in two: update its box and add
                // the new node next to it
                picked.rect = node_cover( *picked.child );
                Branch split_branch;
                split_branch.rect = node_cover( *split_node );
                split_branch.child = std::move( split_node );
                return add_branch( std::move( split_branch ), node );
            }

            /// @return the node created if the given node was split
            std::unique_ptr< Node > add_branch( Branch branch, Node& node )
            {
                if( node.count < max_nodes )
                {
                    node.branches[node.count] = std::move( branch );
                    node.count++;
                    return nullptr;
                }
                return split_node( std::move( branch ), node );
            }

            /// Pick the branch needing the smallest volume increase to contain
            /// the rectangle, the smallest one in case of a tie.
            static index_t pick_branch( const Rect& rect, const Node& node )
            {
                index_t best{ 0 };
                double best_increase{ 0. };
                double best_area{ 0. };
                for( const auto b : Range{ node.count } )
                {
                    const auto& branch_rect = node.branches[b].rect;
                    const auto area = volume( branch_rect );
                    const auto increase =
                        volume( combine( rect, branch_rect ) ) - area;
                    if( b == 0 || increase < best_increase
                        || ( increase == best_increase && area < best_area ) )
                    {
                        best = b;
                        best_area = area;
                        best_increase = increase;
                    }
                }
                return best;
            }

            /// Share the branches of a full node and the new branch between
            /// this node and a new one.
            /// @return the new node
            std::unique_ptr< Node > split_node( Branch branch, Node& node )
            {
                assert( node.count == max_nodes );
                Partition partition;
                for( const auto b : Range{ max_nodes } )
                {
                    partition.branches[b] = std::move( node.branches[b] );
                }
                partition.branches[max_nodes] = std::move( branch );
                choose_partition( partition );

                auto new_node = std::make_unique< Node >();
                new_node->level = node.level;
                node.count = 0;
                std::array< Node*, 2 > target_nodes{ &node, new_node.get() };
                for( const auto b : Range{ max_nodes + 1 } )
                {
                    auto& target = *target_nodes[partition.group[b]];
                    target.branches[target.count] =
                        std::move( partition.branches[b] );
                    target.count++;
                }
                assert( node.count + new_node->count == max_nodes + 1 );
                return new_node;
            }

            /*!
             * Quadratic split: the two branches wasting the most volume
             * together seed the two groups. Then, the remaining branch with
             * the strongest preference for a group joins it, until a group
             * must take all the remaining branches to reach the minimal fill.
             */
            static void choose_partition( Partition& partition )
            {
                constexpr auto total = max_nodes + 1;
                partition.group.fill( NOT_TAKEN );
                pick_seeds( partition );
                while( partition.count[0] + partition.count[1] < total
                       && partition.count[0] < total - MIN_NODES
                       && partition.count[1] < total - MIN_NODES )
                {
                    bool first_time{ true };
                    double biggest_diff{ 0. };
                    index_t chosen{ 0 };
                    index_t better_group{ 0 };
                    for( const auto b : Range{ total } )
                    {
                        if( partition.group[b] != NOT_TAKEN )
                        {
                            continue;
                        }
                        const auto& rect = partition.branches[b].rect;
                        const auto growth0 =
                            volume( combine( rect, partition.cover[0] ) )
                            - partition.area[0];
                        const auto growth1 =
                            volume( combine( rect, partition.cover[1] ) )
                            - partition.area[1];
                        auto diff = growth1 - growth0;
                        index_t group{ 0 };
                        if( diff < 0 )
                        {
                            group = 1;
                            diff = -diff;
                        }
                        if( first_time || diff > biggest_diff )
                        {
                            first_time = false;
                            biggest_diff = diff;
                            chosen = b;
                            better_group = group;
                        }
                        else if( diff == biggest_diff
                                 && partition.count[group]
                                        < partition.count[better_group] )
                        {
                            chosen = b;
                            better_group = group;
                        }
                    }
                    classify( chosen, better_group, partition );
                }
                // One group is full enough: the other takes the rest
                if( partition.count[0] + partition.count[1] < total )
                {
                    const index_t group =
                        partition.count[0] >= total - MIN_NODES ? 1 : 0;
                    for( const auto b : Range{ total } )
                    {
                        if( partition.group[b] == NOT_TAKEN )
                        {
                            classify( b, group, partition );
                        }
                    }
                }
                assert( partition.count[0] >= MIN_NODES
                        && partition.count[1] >= MIN_NODES );
            }

            static void pick_seeds( Partition& partition )
            {
                constexpr auto total = max_nodes + 1;
                std::array< double, total > areas;
                for( const auto b : Range{ total } )
                {
                    areas[b] = volume( partition.branches[b].rect );
                }
                bool first_time{ true };
                double worst{ 0. };
                index_t seed0{ 0 };
                index_t seed1{ 0 };
                for( const auto a : Range{ total - 1 } )
                {
                    for( const auto b : Range{ a + 1, total } )
                    {
                        const auto waste =
                            volume( combine( partition.branches[a].rect,
                                partition.branches[b].rect ) )
                            - areas[a] - areas[b];
                        if( first_time || waste > worst )
                        {
                            first_time = false;
                            worst = waste;
                            seed0 = a;
                            seed1 = b;
                        }
                    }
                }
                classify( seed0, 0, partition );
                classify( seed1, 1, partition );
            }

            static void classify(
                index_t branch, index_t group, Partition& partition )
            {
                assert( partition.group[branch] == NOT_TAKEN );
                partition.group[branch] = group;
                const auto& rect = partition.branches[branch].rect;
                partition.cover[group] =
                    partition.count[group] == 0
                        ? rect
                        : combine( rect, partition.cover[group] );
                partition.area[group] = volume( partition.cover[group] );
                partition.count[group]++;
            }

            /// Replace the branch by the last one of the node
            static void disconnect_branch( Node& node, index_t index )
            {
                assert( index < node.count );
                node.count--;
                if( index != node.count )
                {
                    node.branches[index] =
                        std::move( node.branches[node.count] );
                }
                node.branches[node.count] = Branch{};
            }

            /*!
             * Remove the entry below the node. The underfull nodes on the
             * way back up are detached and stored to reinsert their entries.
             * @return true if the entry was found
             */
            bool remove_rec( const Rect& rect,
                const Data& data,
                Node& node,
                std::vector< std::unique_ptr< Node > >& eliminated_nodes )
            {
                if( node.is_leaf() )
                {
                    for( const auto b : Range{ node.count } )
                    {
                        if( node.branches[b].data == data )
                        {
                            disconnect_branch( node, b );
                            return true;
                        }
                    }
                    return false;
                }
                for( const auto b : Range{ node.count } )
                {
                    auto& branch = node.branches[b];
                    if( !overlap( rect, branch.rect )
                        || !remove_rec(
                            rect, data, *branch.child, eliminated_nodes ) )
                    {
                        continue;
                    }
                    if( branch.child->count >= MIN_NODES )
                    {
                        branch.rect = node_cover( *branch.child );
                    }
                    else
                    {
                        eliminated_nodes.push_back( std::move( branch.child ) );
                        disconnect_branch( node, b );
                    }
                    return true;
                }
                return false;
            }

            /// @return false if the callback stopped the search
            template < typename Callback >
            static bool search_rec(
                const Node& node, const Rect& rect, Callback& callback )
            {
                for( const auto b : Range{ node.count } )
                {
                    const auto& branch = node.branches[b];
                    if( !overlap( rect, branch.rect ) )
                    {
                        continue;
                    }
                    if( node.is_leaf() )
                    {
                        if( !callback( branch.data ) )
                        {
                            return false;
                        }
                    }
                    else if( !search_rec( *branch.child, rect, callback ) )
                    {
                        return false;
                    }
                }
                return true;
            }

        private:
            static constexpr double PI = 3.14159265358979323846;
            static constexpr double UNIT_SPHERE_VOLUME =
                dimension == 2 ? PI : 4. / 3. * PI;

            std::unique_ptr< Node > root_;
        };
    } // namespace internal
} // namespace geode
