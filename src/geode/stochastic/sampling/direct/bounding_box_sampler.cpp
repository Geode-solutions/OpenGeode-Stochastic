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
#include <geode/stochastic/sampling/direct/bounding_box_sampler.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/basic/pimpl_impl.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < index_t dimension >
    class BoundingBoxSampler< dimension >::Impl
    {
    public:
        Impl( const BoundingBox< dimension > &box )
        {
            for( const auto dim : LRange( dimension ) )
            {
                dist_coordinates[dim].min_value = box.min().value( dim );
                dist_coordinates[dim].max_value = box.max().value( dim );
            }
        }
        Point< dimension > sample_uniform( RandomEngine &engine )
        {
            geode::Point< dimension > point;
            for( const auto dim : geode::Range( dimension ) )
            {
                point.set_value(
                    dim, engine.sample_uniform( dist_coordinates[dim] ) );
            }
            return point;
        }

    private:
        std::array< UniformClosed< double >, dimension > dist_coordinates;
    };

    template < index_t dimension >
    BoundingBoxSampler< dimension >::BoundingBoxSampler(
        const BoundingBox< dimension > &box )
        : impl_{ box }
    {
    }

    template < index_t dimension >
    BoundingBoxSampler< dimension >::~BoundingBoxSampler() = default;

    template < index_t dimension >
    Point< dimension > BoundingBoxSampler< dimension >::sample_uniform(
        RandomEngine &engine )
    {
        return impl_->sample_uniform( engine );
    }

    template class opengeode_stochastic_stochastic_api BoundingBoxSampler< 1 >;
    template class opengeode_stochastic_stochastic_api BoundingBoxSampler< 2 >;
    template class opengeode_stochastic_stochastic_api BoundingBoxSampler< 3 >;
} // namespace geode