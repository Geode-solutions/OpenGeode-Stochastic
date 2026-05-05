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
    const Type& ObjectSet< Type >::get_fixed_object( index_t index ) const
    {
        OpenGeodeStochasticStochasticException::check_exception(
            index < fixed_objects_.size(), nullptr,
            OpenGeodeException::TYPE::data,
            "[ObjectSet] - index for fixed object out of range." );
        return fixed_objects_[index];
    }

    template < typename Type >
    const Type& ObjectSet< Type >::get_free_object( index_t index ) const
    {
        OpenGeodeStochasticStochasticException::check_exception(
            index < free_objects_.size(), nullptr,
            OpenGeodeException::TYPE::data,
            "[ObjectSet] - index for free object out of range." );
        return free_objects_[index];
    }

    template < typename Type >
    index_t ObjectSet< Type >::add_fixed_object( Type&& object )
    {
        fixed_objects_.push_back( std::move( object ) );
        return fixed_objects_.size() - 1;
    }

    template < typename Type >
    index_t ObjectSet< Type >::add_free_object( Type&& object )
    {
        free_objects_.push_back( std::move( object ) );
        return free_objects_.size() - 1;
    }

    template < typename Type >
    void ObjectSet< Type >::update_free_object( index_t index, Type&& object )
    {
        OpenGeodeStochasticStochasticException::check_exception(
            index < free_objects_.size(), nullptr,
            OpenGeodeException::TYPE::data,
            "[ObjectSet] - free object index out of range." );
        free_objects_[index] = std::move( object );
    }

    template < typename Type >
    void ObjectSet< Type >::remove_free_object( index_t index )
    {
        const index_t last = free_objects_.size() - 1;
        OpenGeodeStochasticStochasticException::check_exception( index <= last,
            nullptr, OpenGeodeException::TYPE::data,
            "[ObjectSet] - free object index out of range." );
        if( index != last )
        {
            std::swap( free_objects_[index], free_objects_[last] );
        }
        free_objects_.pop_back();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_objects() const
    {
        return free_objects_.size() + fixed_objects_.size();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_fixed_objects() const
    {
        return fixed_objects_.size();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_free_objects() const
    {
        return free_objects_.size();
    }

    template < typename Type >
    bool ObjectSet< Type >::empty() const
    {
        return free_objects_.empty() && fixed_objects_.empty();
    }

    template < typename Type >
    std::string ObjectSet< Type >::string() const
    {
        return absl::StrCat( "ObjectSet ",
            this->name().value_or( this->id().string() ), " (",
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