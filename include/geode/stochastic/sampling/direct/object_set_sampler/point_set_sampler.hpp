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

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/sampling/direct/object_set_sampler/object_set_sampler.hpp>
#include <geode/stochastic/sampling/direct/point_uniform_sampler.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{
    template < index_t dimension >
    class UniformPointSetSampler : public ObjectSetSampler< Point< dimension > >
    {
    public:
        UniformPointSetSampler( const SpatialDomain< dimension >& domain,
            const ObjectSamplerConfig< Point< dimension > >& config )
            : ObjectSetSampler< Point< dimension > >{}, domain_( domain )
        {
            auto volume = domain_.extended_n_volume();
            OpenGeodeStochasticStochasticException::check_exception(
                volume != 0., nullptr, OpenGeodeException::TYPE::data,
                "[UniformPointSetSampler] Undefined Extended Bounding "
                "Box (volume ==0)." );
            this->log_pdf_ = -std::log( volume );
            step_move_ = define_step_for_move( config.move_ratio );
            OpenGeodeStochasticStochasticException::check_exception(
                step_move_ > 0., nullptr, OpenGeodeException::TYPE::data,
                "[UniformPointSetSampler] Undefined step length for move "
                "(value == ",
                step_move_, ")." );
        }

        Point< dimension > sample( RandomEngine& engine ) const override
        {
            return PointUniformSampler::sample< dimension >(
                engine, domain_.extended_box() );
        }

        Point< dimension > change(
            const Point< dimension >& obj, RandomEngine& engine ) const override
        {
            geode::Sphere< dimension > ball{ obj, step_move_ };

            auto new_point =
                PointUniformSampler::sample< dimension >( engine, ball );
            constexpr index_t max_try{ 100 };
            for( const auto try_id : geode::Range{ max_try } )
            {
                geode_unused( try_id );
                if( domain_.extended_contains( new_point ) )
                {
                    return new_point;
                }
                new_point =
                    PointUniformSampler::sample< dimension >( engine, ball );
            }
            throw OpenGeodeStochasticStochasticException{ nullptr,
                OpenGeodeException::TYPE::internal,
                "[UniformPointSetSampler] Cannot find a point in the "
                "extended domain" };
        }

    private:
        double define_step_for_move( double ratio )
        {
            return ratio * domain_.smallest_length();
        }

        bool is_valid_object( const Point< dimension >& obj ) const override
        {
            return domain_.extended_contains( obj );
        }

    private:
        const SpatialDomain< dimension >& domain_;
        double step_move_{ 0. };
    };

    template < index_t dimension >
    struct ObjectSamplerConfig< Point< dimension > >
    {
        // use to define the step for change move (move_ratio*domain volume)
        double move_ratio = 0.1;
    };

    template < index_t dimension >
    std::unique_ptr< ObjectSetSampler< Point< dimension > > >
        build_objectset_sampler( const SpatialDomain< dimension >& domain,
            const ObjectSamplerConfig< Point< dimension > >& config )
    {
        return std::make_unique< UniformPointSetSampler< dimension > >(
            domain, config );
    }

} // namespace geode