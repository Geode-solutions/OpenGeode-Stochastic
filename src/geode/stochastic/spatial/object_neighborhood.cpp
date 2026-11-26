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

#include <geode/stochastic/spatial/object_neighborhood.hpp>

#include <algorithm>
#include <array>

#include <geode/basic/pimpl_impl.hpp>

#include <geode/geometry/bounding_box.hpp>

#include <geode/stochastic/spatial/internal/rtree.hpp>

namespace
{
    template < geode::index_t dimension >
    std::pair< std::array< double, dimension >,
        std::array< double, dimension > >
        bounding_box_bounds( const geode::BoundingBox< dimension >& box )
    {
        std::pair< std::array< double, dimension >,
            std::array< double, dimension > >
            min_max_arrays;
        const auto& min = box.min();
        const auto& max = box.max();
        for( const auto dim_id : geode::LRange{ dimension } )
        {
            min_max_arrays.first[dim_id] = min.value( dim_id );
            min_max_arrays.second[dim_id] = max.value( dim_id );
        }
        return min_max_arrays;
    }
} // namespace

namespace geode
{
    template < index_t dimension >
    class ObjectNeighborhood< dimension >::Impl
    {
    public:
        void add( const BoundingBox< dimension >& box, const ObjectId& id )
        {
            auto [min_bound, max_bound] =
                bounding_box_bounds< dimension >( box );
            tree_.Insert( min_bound.data(), max_bound.data(), id );
        }

        void remove( const BoundingBox< dimension >& box, const ObjectId& id )
        {
            auto [min_bound, max_bound] =
                bounding_box_bounds< dimension >( box );
            tree_.Remove( min_bound.data(), max_bound.data(), id );
        }

        std::vector< ObjectId > get_all_neighbor_ids(
            const BoundingBox< dimension >& box,
            const std::vector< uuid >& targeted_set_ids,
            const std::optional< ObjectId >& exclude_self_id ) const
        {
            auto [min_bound, max_bound] =
                bounding_box_bounds< dimension >( box );
            std::vector< ObjectId > res;
            tree_.Search( min_bound.data(), max_bound.data(),
                [&res, &exclude_self_id, &targeted_set_ids](
                    const ObjectId& cur_id ) -> bool {
                    if( !std::binary_search( targeted_set_ids.begin(),
                            targeted_set_ids.end(), cur_id.set_id ) )
                    {
                        return true;
                    }
                    if( exclude_self_id && exclude_self_id.value() == cur_id )
                    {
                        return true;
                    }
                    res.push_back( cur_id );
                    return true;
                } );
            return res;
        }

    private:
        RTree< ObjectId, double, dimension > tree_;
    };

    template < index_t dimension >
    ObjectNeighborhood< dimension >::ObjectNeighborhood() = default;

    template < index_t dimension >
    ObjectNeighborhood< dimension >::~ObjectNeighborhood() = default;

    template < index_t dimension >
    void ObjectNeighborhood< dimension >::add(
        const BoundingBox< dimension >& box, const ObjectId& id )
    {
        impl_->add( box, id );
    }

    template < index_t dimension >
    void ObjectNeighborhood< dimension >::remove(
        const BoundingBox< dimension >& box, const ObjectId& id )
    {
        impl_->remove( box, id );
    }

    template < index_t dimension >
    void ObjectNeighborhood< dimension >::update(
        const BoundingBox< dimension >& old_box,
        const BoundingBox< dimension >& new_box,
        const ObjectId& id )
    {
        impl_->remove( old_box, id );
        impl_->add( new_box, id );
    }

    template < index_t dimension >
    void ObjectNeighborhood< dimension >::update(
        const BoundingBox< dimension >& box,
        const ObjectId& old_id,
        const ObjectId& new_id )
    {
        impl_->remove( box, old_id );
        impl_->add( box, new_id );
    }

    template < index_t dimension >
    std::vector< ObjectId >
        ObjectNeighborhood< dimension >::get_all_neighbor_ids(
            const BoundingBox< dimension >& box,
            const std::vector< uuid >& targeted_set_ids,
            std::optional< ObjectId > exclude_self_id ) const
    {
        return impl_->get_all_neighbor_ids(
            box, targeted_set_ids, exclude_self_id );
    }

    template class opengeode_stochastic_stochastic_api ObjectNeighborhood< 2 >;
    template class opengeode_stochastic_stochastic_api ObjectNeighborhood< 3 >;

} // namespace geode
