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
#include <geode/stochastic/spatial/details/RTree.hpp>
#include <optional>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( BoundingBox );
}

namespace geode
{
    struct ObjectId
    {
        index_t index;
        uuid set_id;
        bool operator==( const ObjectId& other ) const noexcept
        {
            return index == other.index && set_id == other.set_id;
        }
        bool operator!=( const ObjectId& other ) const noexcept
        {
            return index != other.index || set_id != other.set_id;
        }
    };

    template < index_t dimension >
    class ObjectNeighborhood
    {
    public:
        ObjectNeighborhood();
        ~ObjectNeighborhood() = default;

        void add( const BoundingBox< dimension >& box, const ObjectId& id );
        void update( const BoundingBox< dimension >& old_box,
            const BoundingBox< dimension >& new_box,
            const ObjectId& id );
        void update( const BoundingBox< dimension >& box,
            const ObjectId& old_id,
            const ObjectId& new_id );
        void remove( const BoundingBox< dimension >& box, const ObjectId& id );

        std::vector< ObjectId > get_all_neighbor_ids(
            const BoundingBox< dimension >& box,
            std::optional< ObjectId > exclude_self_id ) const;

    private:
        RTree< ObjectId, double, dimension > tree_;
    };

} // namespace geode