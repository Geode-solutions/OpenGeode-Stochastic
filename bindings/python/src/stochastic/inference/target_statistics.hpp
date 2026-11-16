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

#include <geode/stochastic/inference/target_statistics.hpp>

namespace geode
{
    void define_target_statistics( pybind11::module& module )
    {
        pybind11::class_< TargetStatisticConfig > config( module,
            "TargetStatisticConfig",
            "Expected mean value of a statistic of the model, given by the "
            "name of its term. A simulation is validated when the relative "
            "error |mean - value| / |value| of the tracked mean is below the "
            "tolerance." );
        config.def( pybind11::init( []( std::optional< std::string > term_name,
                                        std::optional< double > value,
                                        std::optional< double > tolerance ) {
            TargetStatisticConfig config{};
            detail::set_if_given( config.term_name, term_name );
            detail::set_if_given( config.value, value );
            detail::set_if_given( config.tolerance, tolerance );
            return config;
        } ),
            pybind11::arg( "term_name" ) = std::nullopt,
            pybind11::arg( "value" ) = std::nullopt, pybind11::kw_only(),
            pybind11::arg( "tolerance" ) = std::nullopt,
            "Create a target (defaults: empty term_name, value 0, "
            "tolerance 0.1)." );
        def_readwrite_named( config, "term_name",
            &TargetStatisticConfig::term_name,
            "Name of the term of the model (see term_names())." );
        def_readwrite_named( config, "value", &TargetStatisticConfig::value,
            "Expected mean of the statistic." );
        def_readwrite_named( config, "tolerance",
            &TargetStatisticConfig::tolerance,
            "Relative tolerance used for validation." );
        config.def( "__repr__", []( const TargetStatisticConfig& self ) {
            return python_call_repr( "TargetStatisticConfig",
                { "term_name=" + python_repr( self.term_name ),
                    "value=" + python_repr( self.value ),
                    "tolerance=" + python_repr( self.tolerance ) } );
        } );
    }
} // namespace geode
