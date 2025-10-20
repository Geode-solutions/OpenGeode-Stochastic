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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace
{
    using namespace geode;

    void test_add_sets_and_objects()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id1 = sets.add_set();
        const auto set_id2 = sets.add_set();

        OPENGEODE_EXCEPTION(
            sets.nb_sets() == 2, "[TestObjectSets] - Expected 2 sets" );

        const auto obj_a =
            sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id1 );
        const auto obj_b =
            sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id1 );
        const auto obj_c =
            sets.add_object( geode::Point2D{ { 5.0, 5.0 } }, set_id2 );

        OPENGEODE_EXCEPTION( sets.nb_objects_in_set( set_id1 ) == 2,
            "[TestObjectSets] - Set 1 should have 2 objects" );
        OPENGEODE_EXCEPTION( sets.nb_objects_in_set( set_id2 ) == 1,
            "[TestObjectSets] - Set 2 should have 1 object" );
        OPENGEODE_EXCEPTION( sets.nb_objects() == 3,
            "[TestObjectSets] - Total object count mismatch" );

        const auto& point = sets.get_object( obj_b );
        const auto result = geode::Point2D{ { 1.0, 1.0 } };
        OPENGEODE_EXCEPTION( point == result,
            "[TestObjectSets] - Wrong object value retrieved" );
    }

    void test_neighbors_by_object_id()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set();

        const auto obj0 =
            sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id );
        const auto obj1 =
            sets.add_object( geode::Point2D{ { 1.0, 0.0 } }, set_id );
        const auto obj2 =
            sets.add_object( geode::Point2D{ { 5.0, 0.0 } }, set_id );

        // Near neighbors should be within distance 2.0
        const auto near_neighbors = sets.neighbors( obj0, 2.0 );

        OPENGEODE_EXCEPTION( !near_neighbors.empty(),
            "[TestObjectSets] - Expected at least one "
            "neighbor within 2.0 distance" );
        OPENGEODE_EXCEPTION( near_neighbors.size() == 1,
            "[TestObjectSets] - Expected exactly one neighbor near obj0" );
        OPENGEODE_EXCEPTION( near_neighbors[0].set_id == obj1.set_id
                                 && near_neighbors[0].index == obj1.index,
            "[TestObjectSets] - Wrong neighbor identified for obj0" );
    }

    void test_neighbors_by_object_value()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set();

        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id );
        sets.add_object( geode::Point2D{ { 1.0, 0.0 } }, set_id );
        sets.add_object( geode::Point2D{ { 3.0, 0.0 } }, set_id );

        const geode::Point2D query{ { 0.5, 0.0 } };
        const auto nearby = sets.neighbors( query, 1.0 );

        OPENGEODE_EXCEPTION( !nearby.empty(),
            "[TestObjectSets] - Expected neighbors near query object" );
        OPENGEODE_EXCEPTION( nearby.size() == 2,
            "[TestObjectSets] - Expected 2 neighbors near query object" );
    }

    void test_string_representation()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set();
        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id );
        sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id );

        const auto desc = sets.string();
        OPENGEODE_EXCEPTION( desc.find( "objects" ) != std::string::npos,
            "[TestObjectSets] - string() output should mention objects" );
    }
} // namespace

int main()
{
    try
    {
        geode::Logger::info( "TEST ObjectSets" );
        geode::StochasticLibrary::initialize();

        test_add_sets_and_objects();
        test_neighbors_by_object_id();
        test_neighbors_by_object_value();
        test_string_representation();

        geode::Logger::info( "TEST ObjectSets SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}