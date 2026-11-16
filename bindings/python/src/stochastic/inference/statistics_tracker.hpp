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
        const auto pyclass_name = absl::StrCat( typestr, "StatisticsTracker" );

        pybind11::class_< Tracker >( module, pyclass_name.c_str(),
            "Mean and variance of the statistics of the model terms over the "
            "realizations of a simulation. It refers to the model of the "
            "runner that created it (it keeps the runner alive)." )
            .def( "statiscal_count", &Tracker::statiscal_count,
                "Return the number of realizations recorded." )

            .def( "means", &Tracker::means, pybind11::return_value_policy::copy,
                "Return the mean of each statistic, in the order of "
                "term_names()." )

            .def( "variances", &Tracker::variances,
                "Return the unbiased variance of each statistic, in the order "
                "of term_names(); 0 when fewer than 2 realizations are "
                "recorded." )

            .def( "term_names", &Tracker::term_names,
                "Return the names of the statistics, in the order of means() "
                "and variances()." )

            .def(
                "means_by_term",
                []( const Tracker& tracker ) {
                    pybind11::dict result;
                    const auto names = tracker.term_names();
                    const auto& means = tracker.means();
                    for( const auto term_id : Indices{ names } )
                    {
                        result[pybind11::str( names[term_id] )] =
                            means[term_id];
                    }
                    return result;
                },
                "Return the mean of each statistic as a dict by term name, in "
                "the order of term_names()." )

            .def(
                "variances_by_term",
                []( const Tracker& tracker ) {
                    pybind11::dict result;
                    const auto names = tracker.term_names();
                    const auto variances = tracker.variances();
                    for( const auto term_id : Indices{ names } )
                    {
                        result[pybind11::str( names[term_id] )] =
                            variances[term_id];
                    }
                    return result;
                },
                "Return the variance of each statistic as a dict by term "
                "name, in the order of term_names()." )

            .def( "__repr__", [pyclass_name]( const Tracker& tracker ) {
                return python_call_repr( pyclass_name,
                    { "statiscal_count="
                            + std::to_string( tracker.statiscal_count() ),
                        "nb_terms="
                            + std::to_string( tracker.means().size() ) } );
            } );
    }

    void define_statistics_tracker( pybind11::module_& module )
    {
        define_statistics_tracker_impl< OwnerSegment2D >( module, "Segment2D" );
    }

} // namespace geode
