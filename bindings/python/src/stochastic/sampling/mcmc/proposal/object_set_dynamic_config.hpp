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

#include "../../../../common.hpp"

#include <geode/stochastic/sampling/mcmc/proposal/object_set_dynamic_config.hpp>

namespace geode
{
    void define_object_set_dynamics( pybind11::module& module )
    {
        pybind11::class_< ObjectSetDynamicsConfig >( module,
            "ObjectSetDynamicsConfig",
            "Relative probabilities of the MCMC moves of an object set." )
            .def( pybind11::init<>() )
            .def_readwrite( "birth_ratio",
                &ObjectSetDynamicsConfig::birth_ratio,
                "Relative probability of birth moves." )
            .def_readwrite( "death_ratio",
                &ObjectSetDynamicsConfig::death_ratio,
                "Relative probability of death moves." )
            .def_readwrite( "change_ratio",
                &ObjectSetDynamicsConfig::change_ratio,
                "Relative probability of change moves." );
    }
} // namespace geode
