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
#include <vector>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/sampling/random_engine.hpp>
#include <geode/stochastic/spatial/object_neighborhood.hpp>

// Random insertions, removals and updates in the neighborhood, each search
// being compared to a brute force search on the same boxes.
namespace
{
    constexpr geode::index_t NB_INITIAL_OBJECTS{ 300 };
    constexpr geode::index_t NB_OPERATIONS{ 2000 };
    constexpr geode::index_t NB_QUERIES_PER_OPERATION{ 3 };
    constexpr double DOMAIN_SIZE{ 100. };
    constexpr double MAX_BOX_SIZE{ 10. };

    template < geode::index_t dimension >
    struct Entry
    {
        geode::BoundingBox< dimension > box;
        geode::ObjectId id;
    };

    double random_value(
        geode::RandomEngine& engine, double min_value, double max_value )
    {
        geode::UniformClosed< double > law;
        law.min_value = min_value;
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
    geode::BoundingBox< dimension > random_box( geode::RandomEngine& engine )
    {
        geode::Point< dimension > min;
        geode::Point< dimension > max;
        // one box in ten is a point, as for point objects
        const auto is_point = random_index( engine, 10 ) == 0;
        for( const auto d : geode::LRange{ dimension } )
        {
            const auto origin = random_value( engine, 0., DOMAIN_SIZE );
            min.set_value( d, origin );
            max.set_value( d,
                is_point ? origin
                         : origin + random_value( engine, 0., MAX_BOX_SIZE ) );
        }
        geode::BoundingBox< dimension > box;
        box.add_point( min );
        box.add_point( max );
        return box;
    }

    template < geode::index_t dimension >
    bool overlap( const geode::BoundingBox< dimension >& box_a,
        const geode::BoundingBox< dimension >& box_b )
    {
        for( const auto d : geode::LRange{ dimension } )
        {
            if( box_a.min().value( d ) > box_b.max().value( d )
                || box_b.min().value( d ) > box_a.max().value( d ) )
            {
                return false;
            }
        }
        return true;
    }

    template < geode::index_t dimension >
    std::vector< geode::ObjectId > brute_force_neighbors(
        const std::vector< Entry< dimension > >& entries,
        const geode::BoundingBox< dimension >& box,
        const std::vector< geode::uuid >& targeted_set_ids,
        const std::optional< geode::ObjectId >& exclude_id )
    {
        std::vector< geode::ObjectId > result;
        for( const auto& entry : entries )
        {
            if( !std::binary_search( targeted_set_ids.begin(),
                    targeted_set_ids.end(), entry.id.set_id ) )
            {
                continue;
            }
            if( exclude_id && exclude_id.value() == entry.id )
            {
                continue;
            }
            if( overlap( box, entry.box ) )
            {
                result.push_back( entry.id );
            }
        }
        std::sort( result.begin(), result.end() );
        return result;
    }

    template < geode::index_t dimension >
    class NeighborhoodChecker
    {
    public:
        explicit NeighborhoodChecker( geode::RandomEngine& engine )
            : engine_( engine ), set_ids_{ geode::uuid{}, geode::uuid{} }
        {
            std::sort( set_ids_.begin(), set_ids_.end() );
        }

        void run()
        {
            for( const auto o : geode::Range{ NB_INITIAL_OBJECTS } )
            {
                geode_unused( o );
                add();
            }
            check_queries();
            for( const auto o : geode::Range{ NB_OPERATIONS } )
            {
                geode_unused( o );
                apply_random_operation();
                check_queries();
            }
            while( !entries_.empty() )
            {
                remove();
            }
            check_queries();
        }

    private:
        void apply_random_operation()
        {
            if( entries_.empty() )
            {
                add();
                return;
            }
            switch( random_index( engine_, 4 ) )
            {
                case 0:
                    add();
                    break;
                case 1:
                    remove();
                    break;
                case 2:
                    move();
                    break;
                default:
                    rename();
                    break;
            }
        }

        void add()
        {
            Entry< dimension > entry{ random_box< dimension >( engine_ ),
                new_id() };
            neighborhood_.add( entry.box, entry.id );
            entries_.push_back( std::move( entry ) );
        }

        void remove()
        {
            const auto position = random_index( engine_, entries_.size() );
            neighborhood_.remove(
                entries_[position].box, entries_[position].id );
            entries_[position] = std::move( entries_.back() );
            entries_.pop_back();
        }

        void move()
        {
            auto& entry = entries_[random_index( engine_, entries_.size() )];
            auto new_box = random_box< dimension >( engine_ );
            neighborhood_.update( entry.box, new_box, entry.id );
            entry.box = std::move( new_box );
        }

        void rename()
        {
            auto& entry = entries_[random_index( engine_, entries_.size() )];
            const auto id = new_id();
            neighborhood_.update( entry.box, entry.id, id );
            entry.id = id;
        }

        geode::ObjectId new_id()
        {
            geode::ObjectId id;
            id.index = next_index_++;
            id.fixed = random_index( engine_, 2 ) == 0;
            id.set_id = set_ids_[random_index( engine_, 2 )];
            return id;
        }

        void check_queries()
        {
            for( const auto q : geode::Range{ NB_QUERIES_PER_OPERATION } )
            {
                geode_unused( q );
                auto box = random_box< dimension >( engine_ );
                box.extends( random_value( engine_, 0., MAX_BOX_SIZE ) );
                std::vector< geode::uuid > targeted_set_ids{ set_ids_[0] };
                if( random_index( engine_, 2 ) == 0 )
                {
                    targeted_set_ids = set_ids_;
                }
                std::optional< geode::ObjectId > exclude_id;
                if( !entries_.empty() && random_index( engine_, 2 ) == 0 )
                {
                    exclude_id =
                        entries_[random_index( engine_, entries_.size() )].id;
                }
                auto result = neighborhood_.get_all_neighbor_ids(
                    box, targeted_set_ids, exclude_id );
                std::sort( result.begin(), result.end() );
                const auto expected = brute_force_neighbors(
                    entries_, box, targeted_set_ids, exclude_id );
                geode::OpenGeodeStochasticStochasticException::test(
                    result == expected, "[TestObjectNeighborhood] ", dimension,
                    "D search found ", result.size(), " objects instead of ",
                    expected.size(), " among ", entries_.size() );
            }
        }

    private:
        geode::RandomEngine& engine_;
        std::vector< geode::uuid > set_ids_;
        geode::ObjectNeighborhood< dimension > neighborhood_;
        std::vector< Entry< dimension > > entries_;
        geode::index_t next_index_{ 0 };
    };

    template < geode::index_t dimension >
    void test_random_operations()
    {
        geode::RandomEngine engine;
        engine.set_seed( 42 );
        NeighborhoodChecker< dimension > checker{ engine };
        checker.run();
    }
} // namespace

int main()
{
    try
    {
        geode::Logger::info( "TEST ObjectNeighborhood" );
        geode::OpenGeodeStochasticStochasticLibrary::initialize();

        test_random_operations< 2 >();
        test_random_operations< 3 >();

        geode::Logger::info( "TEST ObjectNeighborhood SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
