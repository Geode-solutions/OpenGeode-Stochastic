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
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <geode/geometry/point.hpp>

#include <geode/stochastic/sampling/random_engine.hpp>
#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/pairwise_interactions/distance_cutoff.hpp>

// Neighbor search in the contexts of the pairwise terms: every pair of objects
// whose interaction is not null must be found by ObjectSets::neighbors with
// the searching distance of the interaction, for the objects of the state and
// for proposed objects, along random additions, moves and removals.
namespace
{
    constexpr geode::index_t NB_INITIAL_OBJECTS{ 150 };
    constexpr geode::index_t NB_OPERATIONS{ 400 };
    constexpr geode::index_t STATE_CHECK_PERIOD{ 40 };
    constexpr geode::index_t FIXED_OBJECT_PERIOD{ 7 };
    constexpr double DOMAIN_SIZE{ 6. };
    constexpr double MAX_SEGMENT_EXTENT{ 1. };
    constexpr double GRID_STEP{ 0.1 };

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

    template < typename Type, typename Interaction >
    class NeighborSearchChecker
    {
    public:
        NeighborSearchChecker( geode::RandomEngine& engine,
            double cutoff,
            bool on_grid,
            std::string_view name )
            : engine_( engine ),
              interaction_{ cutoff },
              on_grid_{ on_grid },
              name_{ name },
              cutoff_{ cutoff }
        {
            set_ids_.push_back( sets_.add_set( "first" ) );
            set_ids_.push_back( sets_.add_set( "second" ) );
            std::sort( set_ids_.begin(), set_ids_.end() );
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
                check_proposed_object();
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

        void apply_random_operation()
        {
            const auto& set_id = set_ids_[random_index( engine_, 2 )];
            const auto nb_free_objects =
                std::as_const( sets_ ).get_set( set_id ).nb_free_objects();
            const auto operation = random_index( engine_, 3 );
            if( operation == 0 || nb_free_objects == 0 )
            {
                sets_.add_object( random_object(), set_id, false );
                return;
            }
            const geode::ObjectId object_id{
                random_index( engine_, nb_free_objects ), false, set_id
            };
            if( operation == 1 )
            {
                sets_.update_free_object( object_id, random_object() );
                return;
            }
            sets_.remove_free_object( object_id );
        }

        /// Count the interacting objects missed by the neighbors search
        void check_object( const geode::ObjectRef< Type >& object_ref,
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

        /// As PairwiseTerm on the objects of the state
        void check_state()
        {
            for( const auto& object_id : sets_.get_all_object() )
            {
                check_object(
                    { sets_.get_object( object_id ), object_id.set_id },
                    object_id );
            }
        }

        /// As PairwiseTerm on an object proposed by the sampler
        void check_proposed_object()
        {
            const auto object = random_object();
            check_object( { object, set_ids_[0] }, std::nullopt );
        }

    private:
        geode::RandomEngine& engine_;
        Interaction interaction_;
        bool on_grid_;
        std::string_view name_;
        double cutoff_;
        geode::ObjectSets< Type > sets_;
        std::vector< geode::uuid > set_ids_;
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
                NeighborSearchChecker< Type, Interaction > checker{ engine,
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
