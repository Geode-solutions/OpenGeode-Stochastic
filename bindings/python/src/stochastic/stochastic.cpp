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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sampling/direct/double_sampler.hpp"

#include "sampling/mcmc/helpers/fracture_simulation_runner.hpp"
#include "sampling/mcmc/helpers/simulation_monitor.hpp"
#include "sampling/mcmc/helpers/simulation_printer.hpp"
#include "sampling/mcmc/helpers/simulation_runner.hpp"

#include "sampling/distributions.hpp"
#include "sampling/random_engine.hpp"

PYBIND11_MODULE( opengeode_stochastic_py_stochastic, module )
{
    module.doc() = "OpenGeode-Stochastic Python binding";
    pybind11::class_< geode::StochasticLibrary >( module, "StochasticLibrary" )
        .def( "initialize", &geode::StochasticLibrary::initialize );

    geode::define_double_sampler( module );

    geode::define_fracture_simulation( module );
    geode::define_simulation_monitor( module );
    geode::define_simulation_printer( module );
    geode::define_simulation_runner( module );

    geode::define_distributions( module );
    geode::define_random_engine( module );
}