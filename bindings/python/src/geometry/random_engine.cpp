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

#include <geode/stochastic/geometry/random_engine.hpp>

namespace geode
{
    template < typename type >
    type python_sample_uniform( pybind11::class_< RandomEngine >& rand_engine,
        const std::string& suffix )
    {
        const auto read_suffix = absl::StrCat( "sample_uniform_", suffix );
        manager.def(
            read_suffix.c_str(), &RandomEngine::sample_uniform< type > );
    }

    void define_random_engine( pybind11::module& module )
    {
        pybind11::class_< RandomEngine > rand_engine( module, "RandomEngine" );
        rand_engine.def( pybind11::init<>() )
            .def( "set_seed", &RandomEngine::set_seed )
            .def( "sample_gaussian", &RandomEngine::sample_gaussian )
            .def( "sample_bernoulli", &AttributeManager::sample_bernoulli )
            .def( "attribute_exists", &AttributeManager::attribute_exists );
        python_sample_uniform< index_t >( rand_engine, "index_t" );
        python_sample_uniform< local_index_t >( rand_engine, "local_index_t" );
        python_sample_uniform< signed_index_t >(
            rand_engine, "signed_index_t" );
        python_sample_uniform< float >( rand_engine, "float" );
        python_sample_uniform< double >( rand_engine, "double" );
    }
} // namespace geode
