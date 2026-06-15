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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sampling/direct/double_sampler.hpp"

#include "inference/statistics_tracker.hpp"

#include "sampling/mcmc/helpers/simulation_printer.hpp"
#include "sampling/mcmc/simulation_runner.hpp"

#include "sampling/distributions.hpp"
#include "sampling/random_engine.hpp"
#include "spatial/spatial_domain.hpp"

#include "applications/fractures.hpp"

PYBIND11_MODULE( opengeode_stochastic_py_stochastic, module )
{
    module.doc() = "OpenGeode-Stochastic Python binding";
    pybind11::class_< geode::OpenGeodeStochasticStochasticLibrary >(
        module, "OpenGeodeStochasticStochasticLibrary" )
        .def( "initialize",
            &geode::OpenGeodeStochasticStochasticLibrary::initialize );

    geode::define_spatial_domain( module );

    geode::define_distributions( module );
    geode::define_random_engine( module );
    geode::define_double_sampler( module );

    geode::define_statistics_tracker( module );
    geode::define_simulation_printer( module );
    geode::define_simulation_configurator_and_runner( module );

    geode::define_fracture_network_description( module );
}