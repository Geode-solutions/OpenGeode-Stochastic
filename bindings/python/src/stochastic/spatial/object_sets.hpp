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

#include "../../common.hpp"

#include <geode/stochastic/spatial/object_sets.hpp>

namespace geode
{
    template < typename object_type >
    void define_object_sets_impl(
        pybind11::module_& module, const std::string& typestr )
    {
        using ObjectSets = geode::ObjectSets< object_type >;
        const auto pyclass_name = absl::StrCat( typestr, "ObjectSets" );

        pybind11::class_< ObjectSets >( module, pyclass_name.c_str() )
            .def( pybind11::init<>() )
            .def( "nb_sets", &ObjectSets::nb_sets )
            .def( "nb_objects", &ObjectSets::nb_objects );
    }

    void define_object_sets( pybind11::module& module )
    {
        define_object_sets_impl< Point2D >( module, "Point2D" );
        define_object_sets_impl< Point3D >( module, "Point3D" );
        define_object_sets_impl< OwnerSegment2D >( module, "Segment2D" );
    }
} // namespace geode