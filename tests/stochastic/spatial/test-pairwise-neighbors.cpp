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

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <absl/container/flat_hash_map.h>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/models/energy_terms/pairwise_term.hpp>
#include <geode/stochastic/sampling/random_engine.hpp>
#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

// Neighbor search in the contexts of the pairwise terms, compared to a brute
// force over all the pairs of objects, along random additions, moves and
// removals:
// - every pair whose interaction is not null is found by
//   ObjectSets::neighbors, for the objects of the state and for proposed
//   objects;
// - the PairwiseTerm counts each interacting pair once: total_log, statistic
//   and the deltas of the operations match the brute force counts.
// The first set interacts with both sets, the second set does not interact
// with itself.
namespace
{
    constexpr geode::index_t NB_INITIAL_OBJECTS{ 150 };
    constexpr geode::index_t NB_OPERATIONS{ 400 };
    constexpr geode::index_t STATE_CHECK_PERIOD{ 40 };
    constexpr geode::index_t FIXED_OBJECT_PERIOD{ 7 };
    constexpr double DOMAIN_SIZE{ 6. };
    constexpr double VOI_MARGIN{ 1. };
    constexpr double MAX_SEGMENT_EXTENT{ 1. };
    constexpr double GRID_STEP{ 0.1 };
    constexpr double GAMMA{ 0.5 };

    double random_value( geode::RandomEngine& engine, double max_value )
    {
        geode::UniformClosed< double > law;
        law.min_value = 0.;
        law.max_value = max_value;
        return engine.sample_uniform( law );
    }

    geode::index_t random_index(
        geode::RandomEngine& engine, geode::index_t size )
    {
        geode::UniformClosedOpen< double > law;
        law.min_value = 0.;
        law.max_value = static_cast< double >( size );
        return static_cast< geode::index_t >(
            std::floor( engine.sample_uniform( law ) ) );
    }

    template < geode::index_t dimension >
    geode::SpatialDomain< dimension > create_domain()
    {
        geode::Point< dimension > min;
        geode::Point< dimension > max;
        for( const auto d : geode::LRange{ dimension } )
        {
            min.set_value( d, VOI_MARGIN );
            max.set_value( d, DOMAIN_SIZE - VOI_MARGIN );
        }
        geode::BoundingBox< dimension > box;
        box.add_point( min );
        box.add_point( max );
        return geode::SpatialDomain< dimension >{ box, VOI_MARGIN };
    }

    template < typename Type, typename Interaction >
    class PairwiseNeighborsChecker
    {
        using Term = geode::PairwiseTerm< Type >;

    public:
        PairwiseNeighborsChecker( geode::RandomEngine& engine,
            double cutoff,
            bool on_grid,
            std::string_view name )
            : engine_( engine ),
              interaction_{ cutoff },
              on_grid_{ on_grid },
              name_{ name },
              cutoff_{ cutoff },
              domain_{ create_domain< Type::dim >() }
        {
            set_ids_.push_back( sets_.add_set( "first" ) );
            set_ids_.push_back( sets_.add_set( "second" ) );
            auto both_sets = set_ids_;
            std::sort( both_sets.begin(), both_sets.end() );
            adjacency_.emplace( set_ids_[0], both_sets );
            adjacency_.emplace(
                set_ids_[1], std::vector< geode::uuid >{ set_ids_[0] } );
            auto adjacency = adjacency_;
            auto impacted_set_ids = set_ids_;
            term_ = std::make_unique< Term >( "pairwise", GAMMA,
                std::move( impacted_set_ids ), domain_, std::move( adjacency ),
                std::make_unique< Interaction >( cutoff ) );
        }

        void run()
        {
            for( const auto o : geode::Range{ NB_INITIAL_OBJECTS } )
            {
                sets_.add_object( random_object(), set_ids_[o % 2],
                    o % FIXED_OBJECT_PERIOD == 0 );
            }
            check_state();
            for( const auto o : geode::Range{ NB_OPERATIONS } )
            {
                apply_random_operation();
                if( o % STATE_CHECK_PERIOD == 0 )
                {
                    check_state();
                }
            }
            check_state();
            geode::OpenGeodeStochasticStochasticException::test(
                nb_missed_pairs_ == 0, "[TestPairwiseNeighbors] ", name_,
                " with cutoff ", cutoff_, ( on_grid_ ? " on grid" : "" ), ": ",
                nb_missed_pairs_, " missed pairs among ", nb_interacting_pairs_,
                " interacting pairs" );
            geode::OpenGeodeStochasticStochasticException::test(
                nb_interacting_pairs_ > 0, "[TestPairwiseNeighbors] ", name_,
                ": no interacting pair to check" );
        }

    private:
        double random_coordinate( double max_value )
        {
            const auto value = random_value( engine_, max_value );
            // Decimal coordinates give pairs at exactly the cutoff distance
            return on_grid_ ? std::round( value / GRID_STEP ) * GRID_STEP
                            : value;
        }

        template < geode::index_t dimension >
        geode::Point< dimension > random_point()
        {
            geode::Point< dimension > point;
            for( const auto d : geode::LRange{ dimension } )
            {
                point.set_value( d, random_coordinate( DOMAIN_SIZE ) );
            }
            return point;
        }

        Type random_object()
        {
            if constexpr( std::is_same_v< Type, geode::OwnerSegment2D > )
            {
                const auto start = random_point< 2 >();
                auto end = start;
                while( end == start )
                {
                    end = geode::Point2D{ { start.value( 0 )
                                                + random_coordinate(
                                                    2 * MAX_SEGMENT_EXTENT )
                                                - MAX_SEGMENT_EXTENT,
                        start.value( 1 )
                            + random_coordinate( 2 * MAX_SEGMENT_EXTENT )
                            - MAX_SEGMENT_EXTENT } };
                }
                return geode::OwnerSegment2D{ start, end };
            }
            else
            {
                return random_point< Type::dim >();
            }
        }

        /// Check the delta of the operation before applying it
        void apply_random_operation()
        {
            const auto& set_id = set_ids_[random_index( engine_, 2 )];
            const auto nb_free_objects =
                std::as_const( sets_ ).get_set( set_id ).nb_free_objects();
            const auto operation = random_index( engine_, 3 );
            auto object = random_object();
            const geode::ObjectRef< Type > object_ref{ object, set_id };
            check_missed_pairs( object_ref, std::nullopt );
            if( operation == 0 || nb_free_objects == 0 )
            {
                check_value( "delta_log_add",
                    term_->delta_log_add( sets_, object_ref ),
                    term_->contribution(
                        count_interactions( object_ref, std::nullopt ) ) );
                sets_.add_object( std::move( object ), set_id, false );
                return;
            }
            const geode::ObjectId object_id{
                random_index( engine_, nb_free_objects ), false, set_id
            };
            const geode::ObjectRef< Type > old_object_ref{
                sets_.get_object( object_id ), set_id
            };
            const auto old_count =
                count_interactions( old_object_ref, object_id );
            if( operation == 1 )
            {
                const auto new_count =
                    count_interactions( object_ref, object_id );
                check_value( "delta_log_change",
                    term_->delta_log_change( sets_, object_id, object_ref ),
                    term_->contribution( new_count - old_count ) );
                sets_.update_free_object( object_id, std::move( object ) );
                return;
            }
            check_value( "delta_log_remove",
                term_->delta_log_remove( sets_, object_id ),
                term_->contribution( -old_count ) );
            sets_.remove_free_object( object_id );
        }

        [[nodiscard]] bool are_adjacent(
            const geode::uuid& set_a, const geode::uuid& set_b ) const
        {
            const auto& adjacent_sets = adjacency_.at( set_a );
            return std::binary_search(
                adjacent_sets.begin(), adjacent_sets.end(), set_b );
        }

        /// Brute force count of the objects of the adjacent sets interacting
        /// with the given one
        [[nodiscard]] double count_interactions(
            const geode::ObjectRef< Type >& object_ref,
            const std::optional< geode::ObjectId >& excluded_id ) const
        {
            double count{ 0. };
            for( const auto& other_id : sets_.get_all_object() )
            {
                if( ( excluded_id && excluded_id.value() == other_id )
                    || !are_adjacent( object_ref.set_id, other_id.set_id ) )
                {
                    continue;
                }
                count += interaction_.evaluate( object_ref,
                    { sets_.get_object( other_id ), other_id.set_id } );
            }
            return count;
        }

        void check_value(
            std::string_view function, double value, double expected ) const
        {
            geode::OpenGeodeStochasticStochasticException::test(
                value == expected, "[TestPairwiseNeighbors] ", name_,
                " with cutoff ", cutoff_, ( on_grid_ ? " on grid" : "" ), ": ",
                function, " gives ", value, " instead of ", expected );
        }

        /// Count the interacting objects missed by the neighbors search
        void check_missed_pairs( const geode::ObjectRef< Type >& object_ref,
            const std::optional< geode::ObjectId >& object_id )
        {
            auto neighbors = sets_.neighbors( object_ref.object, set_ids_,
                interaction_.neighborhood_searching_distance(), object_id );
            std::sort( neighbors.begin(), neighbors.end() );
            for( const auto& other_id : sets_.get_all_object() )
            {
                if( object_id && object_id.value() == other_id )
                {
                    continue;
                }
                const geode::ObjectRef< Type > other_ref{
                    sets_.get_object( other_id ), other_id.set_id
                };
                if( interaction_.evaluate( object_ref, other_ref ) <= 0. )
                {
                    continue;
                }
                nb_interacting_pairs_++;
                if( !std::binary_search(
                        neighbors.begin(), neighbors.end(), other_id ) )
                {
                    nb_missed_pairs_++;
                }
            }
        }

        /// Neighbors of each object, total_log and statistic of the state
        void check_state()
        {
            const auto object_ids = sets_.get_all_object();
            double nb_pairs{ 0. };
            double nb_pairs_in_domain{ 0. };
            for( const auto a : geode::Indices{ object_ids } )
            {
                const geode::ObjectRef< Type > object_a{
                    sets_.get_object( object_ids[a] ), object_ids[a].set_id
                };
                check_missed_pairs( object_a, object_ids[a] );
                for( const auto b : geode::Range{ a + 1, object_ids.size() } )
                {
                    const geode::ObjectRef< Type > object_b{
                        sets_.get_object( object_ids[b] ), object_ids[b].set_id
                    };
                    if( !are_adjacent( object_a.set_id, object_b.set_id )
                        || interaction_.evaluate( object_a, object_b ) <= 0. )
                    {
                        continue;
                    }
                    nb_pairs++;
                    if( domain_.contains(
                            interaction_.location( object_a, object_b ) ) )
                    {
                        nb_pairs_in_domain++;
                    }
                }
            }
            check_value( "total_log", term_->total_log( sets_ ),
                term_->contribution( nb_pairs ) );
            check_value(
                "statistic", term_->statistic( sets_ ), nb_pairs_in_domain );
        }

    private:
        geode::RandomEngine& engine_;
        Interaction interaction_;
        bool on_grid_;
        std::string_view name_;
        double cutoff_;
        geode::SpatialDomain< Type::dim > domain_;
        geode::ObjectSets< Type > sets_;
        std::vector< geode::uuid > set_ids_;
        absl::flat_hash_map< geode::uuid, std::vector< geode::uuid > >
            adjacency_;
        std::unique_ptr< Term > term_;
        geode::index_t nb_interacting_pairs_{ 0 };
        geode::index_t nb_missed_pairs_{ 0 };
    };

    template < typename Type, typename Interaction >
    void test_interaction( std::string_view name )
    {
        geode::RandomEngine engine;
        engine.set_seed( 42 );
        for( const auto on_grid : { true, false } )
        {
            for( const auto cutoff : { 0.3, 0.7, 1.0 } )
            {
                PairwiseNeighborsChecker< Type, Interaction > checker{ engine,
                    cutoff, on_grid, name };
                checker.run();
            }
        }
    }
} // namespace

int main()
{
    try
    {
        geode::Logger::info( "TEST PairwiseNeighbors" );
        geode::OpenGeodeStochasticStochasticLibrary::initialize();

        test_interaction< geode::Point2D,
            geode::CenterEuclideanDistanceCutoff< geode::Point2D > >(
            "CenterEuclideanDistanceCutoff<Point2D>" );
        test_interaction< geode::Point3D,
            geode::CenterEuclideanDistanceCutoff< geode::Point3D > >(
            "CenterEuclideanDistanceCutoff<Point3D>" );
        test_interaction< geode::OwnerSegment2D,
            geode::CenterEuclideanDistanceCutoff< geode::OwnerSegment2D > >(
            "CenterEuclideanDistanceCutoff<OwnerSegment2D>" );
        test_interaction< geode::Point2D,
            geode::MinimalDistanceCutoff< geode::Point2D > >(
            "MinimalDistanceCutoff<Point2D>" );
        test_interaction< geode::Point3D,
            geode::MinimalDistanceCutoff< geode::Point3D > >(
            "MinimalDistanceCutoff<Point3D>" );
        test_interaction< geode::OwnerSegment2D,
            geode::MinimalDistanceCutoff< geode::OwnerSegment2D > >(
            "MinimalDistanceCutoff<OwnerSegment2D>" );

        geode::Logger::info( "TEST PairwiseNeighbors SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
