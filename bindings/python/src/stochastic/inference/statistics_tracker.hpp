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

#include <geode/stochastic/inference/statistics_tracker.hpp>

#include <geode/geometry/basic_objects/segment.hpp>

namespace geode
{

    template < typename ObjectType >
    void define_statistics_tracker_impl(
        pybind11::module_& module, const std::string& typestr )
    {
        using Tracker = geode::StatisticsTracker< ObjectType >;
        const auto pyclass_name = typestr + "StatisticsTracker";

        pybind11::class_< Tracker >( module, pyclass_name.c_str() )
            .def( "statiscal_count", &Tracker::statiscal_count )

            .def(
                "means", &Tracker::means, pybind11::return_value_policy::copy )

            .def( "variances", &Tracker::variances )

            .def( "__repr__", []( const Tracker& tracker ) {
                return absl::StrCat( "StatisticsTracker(count=",
                    tracker.statiscal_count(), ")" );
            } );
    }

    void define_statistics_tracker( pybind11::module_& module )
    {
        define_statistics_tracker_impl< OwnerSegment2D >(
            module, "OwnerSegment2D" );
    }

} // namespace geode