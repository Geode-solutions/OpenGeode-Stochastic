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
        pybind11::class_< RandomEngine >( module, "RandomEngine" )
            .def( pybind11::init<>() )
            // Seed setting
            .def( "set_seed",
                pybind11::overload_cast< uint64_t >( &RandomEngine::set_seed ),
                pybind11::arg( "number" ), "Set RNG seed using integer" )
            .def( "set_seed",
                pybind11::overload_cast< std::string_view >(
                    &RandomEngine::set_seed ),
                pybind11::arg( "word" ), "Set RNG seed using string" )

            // Uniform sampling (double)
            .def(
                "sample_uniform_closed",
                []( RandomEngine &self, const UniformClosed< double > &law ) {
                    return self.sample_uniform( law );
                },
                pybind11::arg( "law" ),
                "Sample a double from a uniform closed distribution" )

            .def(
                "sample_uniform_closed_open",
                []( RandomEngine &self,
                    const UniformClosedOpen< double > &law ) {
                    return self.sample_uniform( law );
                },
                pybind11::arg( "law" ),
                "Sample a double from a uniform closed-open distribution" )

            // Gaussian
            .def( "sample_gaussian", &RandomEngine::sample_gaussian,
                pybind11::arg( "law" ),
                "Sample a value from a Gaussian distribution" )
            .def( "sample_truncated_gaussian",
                &RandomEngine::sample_truncated_gaussian,
                pybind11::arg( "law" ),
                "Sample a value from a truncated Gaussian" )

            // Other distributions
            .def( "sample_von_mises", &RandomEngine::sample_von_mises,
                pybind11::arg( "law" ),
                "Sample a value from a Von Mises-Fisher" )
            .def( "sample_truncated_lognormal",
                &RandomEngine::sample_truncated_lognormal,
                pybind11::arg( "law" ),
                "Sample a value from a Truncated Log Normal" )
            .def( "sample_truncated_powerlaw",
                &RandomEngine::sample_truncated_powerlaw,
                pybind11::arg( "law" ),
                "Sample a value from a Truncated Power Law" )
            .def( "sample_log", &RandomEngine::sample_log,
                "Return a logarithmically uniform random value" )
            .def( "sample_bernoulli", &RandomEngine::sample_bernoulli,
                pybind11::arg( "probability_of_success" ),
                "Sample a boolean with given success probability" )

            .def( "__repr__", []( const RandomEngine & ) {
                return "<RandomEngine>";
            } );
    }
} // namespace geode
