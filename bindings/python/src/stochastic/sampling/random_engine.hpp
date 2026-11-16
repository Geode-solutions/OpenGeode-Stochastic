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

#include <geode/stochastic/sampling/random_engine.hpp>

namespace geode
{
    void define_random_engine( pybind11::module &module )
    {
        pybind11::class_< RandomEngine >( module, "RandomEngine",
            "Random number generator of the simulations. Every sampling "
            "function raises an OpenGeodeDataError when the parameters of the "
            "law are invalid." )
            .def( pybind11::init<>(),
                "Create a generator with a default (non reproducible) seed; "
                "call set_seed to get a reproducible sequence." )
            .def( "set_seed",
                pybind11::overload_cast< uint64_t >( &RandomEngine::set_seed ),
                pybind11::arg( "number" ),
                "Reseed the generator with an integer: the same number gives "
                "the same sequence of samples." )
            .def( "set_seed",
                pybind11::overload_cast< std::string_view >(
                    &RandomEngine::set_seed ),
                pybind11::arg( "word" ),
                "Reseed the generator from a string. The seed is a hash of "
                "the word that can differ between processes: use an integer "
                "to reproduce a sequence." )
            .def(
                "sample_uniform_closed",
                []( RandomEngine &self, const UniformClosed< double > &law ) {
                    return self.sample_uniform( law );
                },
                pybind11::arg( "law" ),
                "Sample a value of the uniform law on [min_value, max_value]; "
                "requires min_value <= max_value." )
            .def(
                "sample_uniform_closed_open",
                []( RandomEngine &self,
                    const UniformClosedOpen< double > &law ) {
                    return self.sample_uniform( law );
                },
                pybind11::arg( "law" ),
                "Sample a value of the uniform law on [min_value, max_value); "
                "requires min_value < max_value." )
            .def( "sample_gaussian", &RandomEngine::sample_gaussian,
                pybind11::arg( "law" ),
                "Sample a value of the Gaussian law; requires a finite mean "
                "and a finite standard deviation > 0." )
            .def( "sample_truncated_gaussian",
                &RandomEngine::sample_truncated_gaussian,
                pybind11::arg( "law" ),
                "Sample a value of the Gaussian law truncated to [min_value, "
                "max_value] (inverse CDF method); requires min_value < "
                "max_value." )
            .def( "sample_von_mises", &RandomEngine::sample_von_mises,
                pybind11::arg( "law" ),
                "Sample an angle in radians, in [0, 2*pi), of the Von Mises "
                "law." )
            .def( "sample_truncated_lognormal",
                &RandomEngine::sample_truncated_lognormal,
                pybind11::arg( "law" ),
                "Sample a value of the log-normal law truncated to "
                "[min_value, max_value]; requires min_value < max_value." )
            .def( "sample_truncated_powerlaw",
                &RandomEngine::sample_truncated_powerlaw,
                pybind11::arg( "law" ),
                "Sample a value of the power law truncated to [min_value, "
                "max_value]; requires alpha > 0 and min_value < max_value." )
            .def( "sample_log", &RandomEngine::sample_log,
                "Return the natural logarithm of a value drawn uniformly in "
                "(0, 1]: a value <= 0." )
            .def( "sample_bernoulli", &RandomEngine::sample_bernoulli,
                pybind11::arg( "probability_of_success" ),
                "Return True with the given probability, which must be in "
                "[0, 1]." )
            .def( "__repr__", []( const RandomEngine & ) {
                return "RandomEngine()";
            } );
    }
} // namespace geode
