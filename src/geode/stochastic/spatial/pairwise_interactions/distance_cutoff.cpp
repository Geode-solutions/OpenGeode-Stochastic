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

#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/pairwise_interactions.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/distance.hpp>

namespace
{
    template < geode::index_t dimension >
    double compute_center_euclidean_distance(
        const geode::Point< dimension >& point1,
        const geode::Point< dimension >& point2 )
    {
        return geode::point_point_distance( point1, point2 );
    }
    template < geode::index_t dimension >
    double compute_center_euclidean_distance(
        const geode::Segment< dimension >& seg1,
        const geode::Segment< dimension >& seg2 )
    {
        return geode::point_point_distance(
            seg1.barycenter(), seg2.barycenter() );
    }
    template < geode::index_t dimension >
    double compute_min_distance( const geode::Point< dimension >& point1,
        const geode::Point< dimension >& point2 )
    {
        return geode::point_point_distance( point1, point2 );
    }
    template < geode::index_t dimension >
    double compute_min_distance( const geode::Segment< dimension >& seg1,
        const geode::Segment< dimension >& seg2 )
    {
        return std::get< 0 >( geode::segment_segment_distance( seg1, seg2 ) );
    }
} // namespace
namespace geode
{
    template < typename Type >
    CenterEuclideanDistanceCutoff< Type >::CenterEuclideanDistanceCutoff(
        double cutoff_distance )
        : PairwiseInteraction< Type >(), cutoff_distance_( cutoff_distance )
    {
    }

    template < typename Type >
    CenterEuclideanDistanceCutoff< Type >::CenterEuclideanDistanceCutoff(
        double cutoff_distance,
        typename PairwiseInteraction< Type >::SCOPE scope )
        : PairwiseInteraction< Type >( scope ),
          cutoff_distance_( cutoff_distance )
    {
    }

    template < typename Type >
    double
        CenterEuclideanDistanceCutoff< Type >::neighborhood_searching_distance()
            const
    {
        return cutoff_distance_;
    }

    template < typename Type >
    double CenterEuclideanDistanceCutoff< Type >::compute(
        const ObjectRef< Type >& object_a,
        const ObjectRef< Type >& object_b ) const
    {
        auto dist = compute_center_euclidean_distance< Type::dim >(
            object_a.object, object_b.object );
        return dist <= cutoff_distance_ ? 1.0 : 0.0;
    }

    template class opengeode_stochastic_stochastic_api
        CenterEuclideanDistanceCutoff< Point< 2 > >;
    template class opengeode_stochastic_stochastic_api
        CenterEuclideanDistanceCutoff< Point< 3 > >;

    template class opengeode_stochastic_stochastic_api
        CenterEuclideanDistanceCutoff< OwnerSegment< 2 > >;

    template < typename Type >
    MinimalDistanceCutoff< Type >::MinimalDistanceCutoff(
        double cutoff_distance )
        : PairwiseInteraction< Type >(), cutoff_distance_( cutoff_distance )
    {
    }

    template < typename Type >
    MinimalDistanceCutoff< Type >::MinimalDistanceCutoff(
        double cutoff_distance,
        typename PairwiseInteraction< Type >::SCOPE scope )
        : PairwiseInteraction< Type >( scope ),
          cutoff_distance_( cutoff_distance )
    {
    }

    template < typename Type >
    double
        MinimalDistanceCutoff< Type >::neighborhood_searching_distance() const
    {
        return cutoff_distance_;
    }

    template < typename Type >
    double MinimalDistanceCutoff< Type >::compute(
        const ObjectRef< Type >& object_a,
        const ObjectRef< Type >& object_b ) const
    {
        auto dist = compute_min_distance< Type::dim >(
            object_a.object, object_b.object );
        return dist <= cutoff_distance_ ? 1.0 : 0.0;
    }

    template class opengeode_stochastic_stochastic_api
        MinimalDistanceCutoff< Point< 2 > >;
    template class opengeode_stochastic_stochastic_api
        MinimalDistanceCutoff< Point< 3 > >;

    template class opengeode_stochastic_stochastic_api
        MinimalDistanceCutoff< OwnerSegment< 2 > >;
} // namespace geode