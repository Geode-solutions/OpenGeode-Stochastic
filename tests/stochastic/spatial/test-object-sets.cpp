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
#include <geode/geometry/point.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

namespace
{
    using namespace geode;
    void test_add_sets_and_objects()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id1 = sets.add_set( "default_name1" );
        const auto set_id2 = sets.add_set( "default_name2" );

        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_sets() == 2, "[TestObjectSets] - Expected 2 sets" );

        const auto obj_a =
            sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id1, false );
        const auto obj_b =
            sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id1, false );
        const auto obj_c =
            sets.add_object( geode::Point2D{ { 5.0, 5.0 } }, set_id2, false );

        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_objects_in_set( set_id1 ) == 2,
            "[TestObjectSets] - Set 1 should have 2 objects" );
        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_objects_in_set( set_id2 ) == 1,
            "[TestObjectSets] - Set 2 should have 1 object" );
        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_objects() == 3,
            "[TestObjectSets] - Total object count mismatch" );

        const auto& point = sets.get_object( obj_b );
        const auto expected = geode::Point2D{ { 1.0, 1.0 } };
        geode::OpenGeodeStochasticStochasticException::test( point == expected,
            "[TestObjectSets] - Wrong object value retrieved" );
    }
    void test_update_free_object()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set( "default_name" );

        const auto obj =
            sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id, false );

        sets.update_free_object( obj, geode::Point2D{ { 9.0, 9.0 } } );

        const auto& updated = sets.get_object( obj );
        auto result = geode::Point2D{ { 9.0, 9.0 } };
        geode::OpenGeodeStochasticStochasticException::test( updated == result,
            "[TestObjectSets] - Updating free object failed" );
    }
    void test_remove_objects()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set( "default_name" );

        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id, false );
        sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id, false );
        sets.add_object( geode::Point2D{ { 2.0, 2.0 } }, set_id, false );

        sets.remove_free_object( { 1, false, set_id } );

        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_objects_in_set( set_id ) == 2,
            "[TestObjectSets] - Expected 2 objects after free removal" );

        // Now remove the last remaining free object using remove_object()
        sets.remove_free_object( { 1, false, set_id } );

        geode::OpenGeodeStochasticStochasticException::test(
            sets.nb_objects_in_set( set_id ) == 1,
            "[TestObjectSets] - Expected 1 object after second removal" );
    }
    void test_get_all_objects()
    {
        ObjectSets< geode::Point2D > sets;
        const auto idA = sets.add_set( "A" );
        const auto idB = sets.add_set( "B" );

        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, idA, false );
        sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, idB, false );

        const auto all = sets.get_all_object();

        geode::OpenGeodeStochasticStochasticException::test( all.size() == 2,
            "[TestObjectSets] - get_all_objects size mismatch" );
    }

    void test_neighbors_by_object_id()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set( "default_name" );

        const auto obj0 =
            sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id, false );
        const auto obj1 =
            sets.add_object( geode::Point2D{ { 1.0, 0.0 } }, set_id, false );
        const auto obj2 =
            sets.add_object( geode::Point2D{ { 5.0, 0.0 } }, set_id, false );

        std::vector< uuid > targeted_set_ids{ set_id };

        const auto near_neighbors = sets.neighbors(
            sets.get_object( obj0 ), targeted_set_ids, 2.0, obj0 );

        geode::OpenGeodeStochasticStochasticException::test(
            near_neighbors.size() == 1,
            "[TestObjectSets] - Expected exactly one neighbor near obj0" );
        geode::OpenGeodeStochasticStochasticException::test(
            near_neighbors[0].set_id == obj1.set_id
                && near_neighbors[0].index == obj1.index,
            "[TestObjectSets] - Wrong neighbor identified for obj0" );
    }

    void test_neighbors_by_object_value()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set( "default_name" );

        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id, false );
        sets.add_object( geode::Point2D{ { 1.0, 0.0 } }, set_id, false );
        sets.add_object( geode::Point2D{ { 3.0, 0.0 } }, set_id, true );

        const geode::Point2D query{ { 0.5, 0.0 } };
        std::vector< uuid > targeted_set_ids{ set_id };
        const auto nearby =
            sets.neighbors( query, targeted_set_ids, 1.0, std::nullopt );

        geode::OpenGeodeStochasticStochasticException::test( nearby.size() == 2,
            "[TestObjectSets] - Expected 2 neighbors near query object" );
    }

    void test_string_representation()
    {
        ObjectSets< geode::Point2D > sets;
        const auto set_id = sets.add_set( "default_name" );
        sets.add_object( geode::Point2D{ { 0.0, 0.0 } }, set_id, false );
        sets.add_object( geode::Point2D{ { 1.0, 1.0 } }, set_id, false );

        const auto desc = sets.string();
        geode::OpenGeodeStochasticStochasticException::test(
            desc.find( "objects" ) != std::string::npos,
            "[TestObjectSets] - string() output should mention objects" );
    }
} // namespace

int main()
{
    try
    {
        geode::Logger::info( "TEST ObjectSets" );
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        test_add_sets_and_objects();
        test_update_free_object();
        test_remove_objects();
        test_get_all_objects();
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
