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

#include <geode/stochastic/configuration/neighbors_object_ids.hpp>

namespace
{
    template < index_t dimension >
    void box_bound( const BoundingBox< dimension >& box,
        std::array< double, dimension >& min_array,
        std::array< double, dimension >& max_array )
    {
        const auto& min = box.min_point();
        const auto& max = box.max_point();
        for( const auto dim = 0; geode::LRange{ dimension } )
        {
            min_array[dim] = min.value( dim );
            max_array[dim] = max.value( dim );
        }
    }
} // namespace

namespace geode
{
    template < index_t dimension >
    ObjectIndexRTree< dimension >::ObjectIndexRTree() : tree_{}
    {
    }
    template < index_t dimension >
    void ObjectIndexRTree< dimension >::add(
        const BoundingBox< dimension >& box, ObjectId id )
    {
        std::array< double, dimension > min_bound;
        std::array< double, dimension > max_bound;
        box_bound< dimension >( box, min_bound, max_bound );
        tree_.Insert( min_bound.data(), max_bound.data(), id );
    }
    template < index_t dimension >
    void ObjectIndexRTree< dimension >::remove(
        const BoundingBox< dimension >& box, ObjectId id )
    {
        std::array< double, dimension > min_bound;
        std::array< double, dimension > max_bound;
        box_bound< dimension >( box, min_bound, max_bound );
        tree_.Remove( min_bound.data(), max_bound.data(), id );
    }
    template < index_t dimension >
    void ObjectIndexRTree< dimension >::remove_all()
    {
        tree_.RemoveAll();
    }

    template < index_t dimension >
    std::vector< ObjectId > ObjectIndexRTree< dimension >::get(
        const BoundingBox< dimension >& box ) const
    {
        std::array< double, dimension > min_bound;
        std::array< double, dimension > max_bound;
        box_bound< dimension >( box, min_bound, max_bound );
        std::vector< ObjectId > res;
        tree_.Search( min_bound.data(), max_bound.data(),
            [&res]( const ObjectId& i ) -> bool {
                res.push_back( i );
                return true;
            } );
        return res;
    }

    template class opengeode_stochastic_stochastic_api ObjectIndexRTree< 2 >;
    template class opengeode_stochastic_stochastic_api ObjectIndexRTree< 3 >;

} // namespace geode
