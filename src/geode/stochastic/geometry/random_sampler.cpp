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

#include <geode/stochastic/geometry/random_sampler.hpp>

#include <geode/stochastic/geometry/common.hpp>
#include <geode/stochastic/geometry/random_engine.hpp>

#include <geode/geometry/point.hpp>

namespace
{
    template < geode::index_t dimension >
    geode::Point< dimension > sample_point_in_box( geode::RandomEngine& engine,
        const geode::UniformBox< dimension >& dist )
    {
        geode::Point< dimension > point;
        for( const auto i : geode::Range( dimension ) )
        {
            point.set_value( i, geode::RandomSampler< double >::sample(
                                    engine, dist.dist_coordinates[i] ) );
        }
        return point;
    }
    geode::Point2D sample_point_in_ball(
        geode::RandomEngine& engine, const geode::UniformBall< 2 >& ball )
    {
        geode::Point< 2 > point;
        geode::UniformClosed< double > dist;
        dist.min_value = 0.;
        dist.max_value = 1.;
        double r =
            ball.radius * std::sqrt( engine.sample_uniform< double >( dist ) );
        double theta = ( 2.0 * M_PI ) * engine.sample_uniform< double >( dist );
        point.set_value( 0, ball.center.value( 0 ) + r * cos( theta ) );
        point.set_value( 1, ball.center.value( 1 ) + r * sin( theta ) );
        return point;
    }
    geode::Point3D sample_point_in_ball(
        geode::RandomEngine& engine, const geode::UniformBall< 3 >& ball )
    {
        geode::Point< 3 > point;
        geode::UniformClosed< double > dist;
        dist.min_value = 0.;
        dist.max_value = 1.;
        auto u = 2.0 * engine.sample_uniform< double >( dist ) - 1.0;
        auto v = std::sqrt( 1. - std::pow( u, 2.0 ) );
        auto theta = 2.0 * M_PI * engine.sample_uniform< double >( dist );

        point.set_value(
            0, ball.center.value( 0 ) + ball.radius * v * cos( theta ) );
        point.set_value(
            1, ball.center.value( 1 ) + ball.radius * v * sin( theta ) );
        point.set_value( 2, ball.center.value( 2 ) + ball.radius * u );

        return point;
    }
} // namespace
namespace geode
{
    double RandomSampler< double >::sample(
        RandomEngine& engine, const DoubleDistribution& dist )
    {
        return std::visit(
            [&engine]( auto&& d ) {
                using D = std::decay_t< decltype( d ) >;
                if constexpr( std::is_same_v< D, UniformClosed< double > > )
                    return engine.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, UniformClosedOpen< double > > )
                    return engine.sample_uniform< double >( d );
                if constexpr( std::is_same_v< D, Gaussian > )
                    return engine.sample_gaussian( d );
                if constexpr( std::is_same_v< D, TruncatedGaussian > )
                    return engine.sample_truncated_gaussian( d );
                OPENGEODE_EXCEPTION(
                    "RandomSampler - Unsupported distribution for double" );
            },
            dist );
    }

    Point2D RandomSampler< Point2D >::sample(
        RandomEngine& engine, const SpatialDistribution2D& dist )
    {
        return std::visit(
            [&engine]( auto&& d ) {
                using D = std::decay_t< decltype( d ) >;
                if constexpr( std::is_same_v< D, UniformBox< 2 > > )
                    return sample_point_in_box( engine, d );
                if constexpr( std::is_same_v< D, UniformBall< 2 > > )
                    return sample_point_in_ball( engine, d );
                OPENGEODE_EXCEPTION( "RandomSampler - Unsupported spatial "
                                     "distribution for Point2D" );
            },
            dist );
    }
    Point3D RandomSampler< Point3D >::sample(
        RandomEngine& engine, const SpatialDistribution3D& dist )
    {
        return std::visit(
            [&engine]( auto&& d ) {
                using D = std::decay_t< decltype( d ) >;
                if constexpr( std::is_same_v< D, geode::UniformBox< 3 > > )
                    return sample_point_in_box( engine, d );
                if constexpr( std::is_same_v< D, geode::UniformBall< 3 > > )
                    return sample_point_in_ball( engine, d );
                OPENGEODE_EXCEPTION( "RandomSampler - Unsupported spatial "
                                     "distribution for Point3D" );
            },
            dist );
    }

} // namespace geode