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
#include <geode/stochastic/sampling/direct/ball_sampler.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/basic/pimpl_impl.hpp>

#include <geode/geometry/basic_objects/sphere.hpp>
#include <geode/geometry/point.hpp>

namespace
{
    geode::Point2D sample_point_in_ball( geode::RandomEngine& engine,
        const geode::UniformClosed< double >& uniform_dist,
        const geode::Point2D& center,
        double radius )
    {
        geode::Point2D point;
        double r =
            radius
            * std::sqrt( engine.sample_uniform< double >( uniform_dist ) );
        double theta =
            ( 2.0 * M_PI ) * engine.sample_uniform< double >( uniform_dist );
        point.set_value( 0, center.value( 0 ) + r * cos( theta ) );
        point.set_value( 1, center.value( 1 ) + r * sin( theta ) );
        return point;
    }

    geode::Point3D sample_point_in_ball( geode::RandomEngine& engine,
        const geode::UniformClosed< double >& uniform_dist,
        const geode::Point3D& center,
        double radius )
    {
        geode::Point3D point;
        auto u = 2.0 * engine.sample_uniform< double >( uniform_dist ) - 1.0;
        auto v = std::sqrt( 1. - std::pow( u, 2.0 ) );
        auto theta =
            2.0 * M_PI * engine.sample_uniform< double >( uniform_dist );

        point.set_value( 0, center.value( 0 ) + radius * v * cos( theta ) );
        point.set_value( 1, center.value( 1 ) + radius * v * sin( theta ) );
        point.set_value( 2, center.value( 2 ) + radius * u );

        return point;
    }
} // namespace
namespace geode
{

    template < index_t dimension >
    class BallSampler< dimension >::Impl
    {
    public:
        Impl( const Sphere< dimension >& ball ) : ball_{ ball }
        {
            dist_.min_value = 0.;
            dist_.max_value = 1.;
        }

        geode::Point< dimension > sample_uniform( geode::RandomEngine& engine )
        {
            return sample_point_in_ball(
                engine, dist_, ball_.origin(), ball_.radius() );
        }

    private:
        const Sphere< dimension >& ball_;
        geode::UniformClosed< double > dist_;
    };

    template < index_t dimension >
    BallSampler< dimension >::BallSampler( const Sphere< dimension >& ball )
        : impl_( ball ){};

    template < index_t dimension >
    BallSampler< dimension >::~BallSampler() = default;

    template < index_t dimension >
    Point< dimension > BallSampler< dimension >::sample_uniform(
        RandomEngine& engine )
    {
        return impl_->sample_uniform( engine );
    }

    template class opengeode_stochastic_stochastic_api BallSampler< 2 >;
    template class opengeode_stochastic_stochastic_api BallSampler< 3 >;
} // namespace geode