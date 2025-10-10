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

#pragma once

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>

namespace geode
{
    template < index_t dimension >
    class UniformPointSetSampler : public ObjectSetSampler< Point< dimension > >
    {
    public:
        UniformPointSetSampler( const BoundingBox< dimension >& box )
            : ObjectSetSampler< Point< dimension > >{}, box_( box )
        {
            auto volume = box_.n_volume();
            OPENGEODE_EXCEPTION( volume != 0.,
                "[PointSetSampler] - Undefined Bounding Box (volume ==0)." );
            this->log_pdf_ = -std::log( volume );
        }

        Point< dimension > sample( RandomEngine& engine ) const override
        {
            return PointUniformSampler::sample< dimension >( engine, box_ );
        }

        Point< dimension > change(
            const Point< dimension >& obj, RandomEngine& engine ) const override
        {
            double ratio = 0.1;
            geode::Sphere< dimension > ball{ obj,
                ratio * std::get< 1 >( box_.smallest_length() ) };

            auto new_point =
                PointUniformSampler::sample< dimension >( engine, ball );
            constexpr index_t max_try{ 100 };
            for( const auto try_id : geode::Range{ max_try } )
            {
                if( box_.contains( new_point ) )
                {
                    return new_point;
                }
                new_point =
                    PointUniformSampler::sample< dimension >( engine, ball );
            }
            throw OpenGeodeException( absl::StrCat(
                "[PointSampler] - Cannot find a point in the box: ",
                box_.string() ) );
            return obj;
        }

    private:
        bool is_valid_object( const Point< dimension >& obj ) const override
        {
            return box_.contains( obj );
        }

    private:
        const BoundingBox< dimension >& box_;
    };

} // namespace geode