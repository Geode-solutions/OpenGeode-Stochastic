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

#include <geode/stochastic/sampling/direct/configuration_sampler/configuration_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>

namespace geode
{
    template < index_t dimension >
    class UniformPointConfigurationSampler
        : public ConfigurationSampler< Point< dimension > >
    {
    public:
        UniformPointConfigurationSampler(
            const BoundingBox< dimension >& box, uuid group_id )
            : ConfigurationSampler< Point< dimension > >{ group_id },
              box_( box )
        {
            auto volume = box_.n_volume();
            if( volume != 0. )
            {
                log_pdf_ = -std::log( volume );
            }
        }

        std::pair< Point< dimension >, uuid > sample(
            RandomEngine& engine ) const override
        {
            return { PointUniformSampler::sample< dimension >( engine, box_ ),
                this->group_id_ };
        }

        std::pair< Point< dimension >, uuid > change(
            const Point< dimension >& obj, RandomEngine& engine ) const override
        {
            double ratio = 0.1;
            geode::Sphere< dimension > ball{ obj,
                ratio * std::get< 1 >( box_.smallest_length() ) };

            auto new_point =
                PointUniformSampler::sample< dimension >( engine, ball );
            while( !box_.contains( new_point ) )
            {
                new_point =
                    PointUniformSampler::sample< dimension >( engine, ball );
            }
            return { new_point, this->group_id_ };
        }

        double log_pdf( const Point< dimension >& obj ) const override
        {
            if( !box_.contains( obj ) )
            {
                return -std::numeric_limits< double >::infinity();
            }
            return log_pdf_;
        }

    private:
        const BoundingBox< dimension >& box_;
        double log_pdf_{ -std::numeric_limits< double >::infinity() };
    };

} // namespace geode