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

#include <cmath>
#include <optional>
#include <string>
#include <string_view>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <absl/container/fixed_array.h>
#include <absl/container/flat_hash_map.h>
#include <absl/container/inlined_vector.h>
#include <absl/types/span.h>

namespace pybind11::detail
{

    template < typename Type >
    struct type_caster< absl::FixedArray< Type > >
        : list_caster< absl::FixedArray< Type >, Type >
    {
    };

    template < typename Type, size_t dimension >
    struct type_caster< absl::InlinedVector< Type, dimension > >
        : list_caster< absl::InlinedVector< Type, dimension >, Type >
    {
    };

    template < typename Type >
    struct type_caster< absl::Span< Type > >
        : list_caster< absl::Span< Type >, Type >
    {
        using value_conv = make_caster< Type >;

        bool load( handle src, bool convert )
        {
            cpp_.clear();
            auto s = reinterpret_borrow< sequence >( src );
            cpp_.reserve( s.size() );
            for( auto it : s )
            {
                value_conv conv;
                if( !conv.load( it, convert ) )
                    return false;
                cpp_.push_back( cast_op< Type&& >( std::move( conv ) ) );
            }
            this->value = absl::MakeConstSpan( cpp_ );
            return true;
        }

        std::vector< typename std::remove_const< Type >::type > cpp_;
    };

    template < typename Key, typename Value >
    struct type_caster< absl::flat_hash_map< Key, Value > >
        : map_caster< absl::flat_hash_map< Key, Value >, Key, Value >
    {
    };
} // namespace pybind11::detail

namespace geode
{
    /*!
     * Return the Python expression of a double: its repr, except for the
     * values that have no literal, which are given by float(...)
     */
    inline std::string python_repr( double value )
    {
        if( std::isnan( value ) )
        {
            return "float('nan')";
        }
        if( std::isinf( value ) )
        {
            return value > 0 ? "float('inf')" : "float('-inf')";
        }
        return pybind11::repr( pybind11::float_( value ) )
            .template cast< std::string >();
    }

    /*!
     * Same as above for an optional double: None when it has no value.
     */
    inline std::string python_repr( const std::optional< double >& value )
    {
        return value ? python_repr( *value ) : "None";
    }

    /*!
     * Return the Python repr of the value, e.g. 1.0 for a double, 'text' for
     * a string or Point2D([1.0, 2.0]) for a Point2D.
     */
    template < typename Type >
    std::string python_repr( const Type& value )
    {
        return pybind11::repr( pybind11::cast( value ) )
            .template cast< std::string >();
    }

    /*!
     * Return name(arguments...), where each argument is already a Python
     * expression (a positional value or keyword=value).
     */
    inline std::string python_call_repr(
        std::string_view name, std::initializer_list< std::string > arguments )
    {
        std::string result{ name };
        result += '(';
        auto first = true;
        for( const auto& argument : arguments )
        {
            if( !first )
            {
                result += ", ";
            }
            first = false;
            result += argument;
        }
        result += ')';
        return result;
    }

    /*!
     * Same as pybind11::class_::def_readwrite, but the setter argument is
     * called "value" instead of "arg0" in the Python signatures.
     */
    template < typename Class,
        typename Member,
        typename Owner,
        typename... Options >
    void def_readwrite_named(
        pybind11::class_< Class, Options... >& python_class,
        const char* name,
        Member Owner::* member,
        const char* doc )
    {
        const pybind11::cpp_function getter(
            [member]( const Class& object ) -> const Member& {
                return object.*member;
            },
            pybind11::is_method( python_class ) );
        const pybind11::cpp_function setter(
            [member]( Class& object, const Member& value ) {
                object.*member = value;
            },
            pybind11::is_method( python_class ), pybind11::arg( "value" ) );
        python_class.def_property( name, getter, setter,
            pybind11::return_value_policy::reference_internal, doc );
    }
} // namespace geode

namespace geode::detail
{
    /// Used by the keyword constructors: fields keep their C++ default value
    /// when the argument is not given.
    template < typename Field, typename Value >
    void set_if_given( Field& field, const std::optional< Value >& value )
    {
        if( value )
        {
            field = *value;
        }
    }
} // namespace geode::detail
