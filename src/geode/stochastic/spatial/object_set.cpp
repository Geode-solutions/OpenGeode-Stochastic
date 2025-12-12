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

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

#include <geode/basic/identifier_builder.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/stochastic/spatial/object_set.hpp>

namespace geode
{
    template < typename Type >
    void ObjectSet< Type >::set_name( std::string_view name )
    {
        IdentifierBuilder builder( *this );
        builder.set_name( name );
    }
    template < typename Type >
    const Type& ObjectSet< Type >::get_object( index_t index ) const
    {
        OPENGEODE_EXCEPTION( index < objects_.size(),
            "[ObjectSet] - object index out of range." );
        return objects_[index];
    }

    template < typename Type >
    index_t ObjectSet< Type >::add_fixed_object( Type&& object )
    {
        auto new_fixed_object_id = first_free_object_;
        first_free_object_++;
        objects_.push_back( std::move( object ) );
        if( new_fixed_object_id != objects_.size() - 1 )
        {
            std::swap( objects_[new_fixed_object_id], objects_.back() );
        }
        return new_fixed_object_id;
    }

    template < typename Type >
    index_t ObjectSet< Type >::add_free_object( Type&& object )
    {
        objects_.push_back( std::move( object ) );
        return objects_.size() - 1;
    }

    template < typename Type >
    void ObjectSet< Type >::update_free_object( index_t index, Type&& object )
    {
        OPENGEODE_EXCEPTION( index >= first_free_object_,
            "[ObjectSet] - cannot update fixed object." );
        OPENGEODE_EXCEPTION( index < objects_.size(),
            "[ObjectSet] - object index out of range." );
        objects_[index] = std::move( object );
    }

    template < typename Type >
    void ObjectSet< Type >::remove_free_object( index_t index )
    {
        OPENGEODE_EXCEPTION( index >= first_free_object_,
            "[ObjectSet] - cannot remove fixed object." );
        remove_object( index );
    }

    template < typename Type >
    void ObjectSet< Type >::remove_object( index_t index )
    {
        const index_t last = objects_.size() - 1;
        OPENGEODE_EXCEPTION(
            index <= last, "[ObjectSet] - object index out of range." );
        if( index != last )
        {
            std::swap( objects_[index], objects_[last] );
        }
        objects_.pop_back();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_objects() const
    {
        return objects_.size();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_fixed_objects() const
    {
        return first_free_object_;
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_free_objects() const
    {
        return objects_.size() - first_free_object_;
    }

    template < typename Type >
    bool ObjectSet< Type >::empty() const
    {
        return objects_.empty();
    }

    template < typename Type >
    bool ObjectSet< Type >::is_fixed( index_t index ) const
    {
        return index < first_free_object_;
    }

    template < typename Type >
    std::string ObjectSet< Type >::string() const
    {
        return absl::StrCat( "ObjectSet ", this->name(), " (",
            this->id().string(), ") contains ", nb_objects(),
            " objects (fixed: ", nb_fixed_objects(),
            " - free: ", nb_free_objects(), ")" );
    }

    // Explicit template instantiation
    template class opengeode_stochastic_stochastic_api ObjectSet< Point2D >;
    template class opengeode_stochastic_stochastic_api ObjectSet< Point3D >;
    template class opengeode_stochastic_stochastic_api
        ObjectSet< OwnerSegment2D >;
} // namespace geode