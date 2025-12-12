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
#include <geode/stochastic/spatial/object_set.hpp>

namespace
{
    void test_add_and_access()
    {
        geode::ObjectSet< geode::Point2D > set;

        OPENGEODE_EXCEPTION(
            set.empty(), "[TestObjectSet] - Set should start empty" );

        const auto idx0 =
            set.add_fixed_object( geode::Point2D{ { 0.0, 0.0 } } );
        const auto idx1 =
            set.add_fixed_object( geode::Point2D{ { 1.0, 1.0 } } );
        const auto idx2 =
            set.add_fixed_object( geode::Point2D{ { 2.0, 2.0 } } );

        OPENGEODE_EXCEPTION( set.nb_objects() == 3,
            "[TestObjectSet] - Set size should be 3 after insertions" );
        OPENGEODE_EXCEPTION( !set.empty(),
            "[TestObjectSet] - Set should not be empty after insertions" );

        const auto& p = set.get_object( idx1 );
        const auto result = geode::Point2D{ { 1.0, 1.0 } };
        OPENGEODE_EXCEPTION(
            p == result, "[TestObjectSet] - Wrong object value at index 1" );
    }

    void test_update_object()
    {
        geode::ObjectSet< geode::Point2D > set;
        const auto idx0 =
            set.add_fixed_object( geode::Point2D{ { 0.0, 0.0 } } );
        const auto idx1 = set.add_free_object( geode::Point2D{ { 1.0, 1.0 } } );

        set.update_free_object( idx1, geode::Point2D{ { 5.0, 5.0 } } );

        const auto& updated = set.get_object( idx1 );
        const auto new_point = geode::Point2D{ { 5.0, 5.0 } };
        OPENGEODE_EXCEPTION(
            updated == new_point, "[TestObjectSet] - Object update failed" );
    }

    void test_remove_object()
    {
        geode::ObjectSet< geode::Point2D > set;
        set.add_fixed_object( geode::Point2D{ { 0.0, 0.0 } } );
        set.add_free_object( geode::Point2D{ { 1.0, 1.0 } } );
        set.add_free_object( geode::Point2D{ { 2.0, 2.0 } } );

        set.remove_free_object( 1 ); // remove the second object (free)

        OPENGEODE_EXCEPTION( set.nb_objects() == 2,
            "[TestObjectSet] - Set size should be 2 after removal" );

        const auto& last = set.get_object( 1 );
        const auto result = geode::Point2D{ { 2.0, 2.0 } };
        OPENGEODE_EXCEPTION( last == result,
            "[TestObjectSet] - Remaining objects not shifted "
            "properly after removal" );
    }

    void test_const_access()
    {
        geode::ObjectSet< geode::Point2D > set;
        set.add_fixed_object( geode::Point2D{ { 10.0, 10.0 } } );

        const auto& const_set = set;
        const auto& p = const_set.get_object( 0 );
        const auto result = geode::Point2D{ { 10.0, 10.0 } };
        OPENGEODE_EXCEPTION(
            p == result, "[TestObjectSet] - Const access mismatch" );
    }

    void test_string_representation()
    {
        geode::ObjectSet< geode::Point2D > set;
        set.add_fixed_object( geode::Point2D{ { 0.0, 0.0 } } );
        set.add_free_object( geode::Point2D{ { 1.0, 1.0 } } );

        const auto desc = set.string();
        OPENGEODE_EXCEPTION( desc.find( "2 objects" ) != std::string::npos,
            "[TestObjectSet] - string() output incorrect" );
    }

    void test_fixed_and_free_management()
    {
        geode::ObjectSet< geode::Point2D > set;

        // Add 2 fixed objects
        auto f0 = set.add_fixed_object( geode::Point2D{ { 0., 0. } } );
        auto u0 = set.add_free_object( geode::Point2D{ { 2., 2. } } );
        auto f1 = set.add_fixed_object( geode::Point2D{ { 1., 1. } } );
        // ici u0 est changé? doije retourner l indice?

        OPENGEODE_EXCEPTION( set.nb_fixed_objects() == 2, "2 fixed objects" );
        OPENGEODE_EXCEPTION( set.nb_free_objects() == 1, "one free object" );
        OPENGEODE_EXCEPTION( f0 == 0 && f1 == 1 );

        // Add 2 free objects
        auto u1 = set.add_free_object( geode::Point2D{ { 2., 2. } } );
        auto u2 = set.add_free_object( geode::Point2D{ { 3., 3. } } );

        OPENGEODE_EXCEPTION(
            set.nb_fixed_objects() == 2, "still two fixed object" );
        OPENGEODE_EXCEPTION( set.nb_free_objects() == 3, "three free objects" );
        OPENGEODE_EXCEPTION( u1 == 3 && u2 == 4 );

        // Remove first fixed object
        set.remove_fixed_object( 0 );

        // Invariant
        OPENGEODE_EXCEPTION(
            set.nb_fixed_objects() == 1, "one remaining fixed objects" );
        OPENGEODE_EXCEPTION(
            set.nb_free_objects() == 3, "still 3 free objects" );
        OPENGEODE_EXCEPTION( set.nb_objects() == 4, "4 total objects" );

        // Remove a free object
        set.remove_free_object( 2 );

        OPENGEODE_EXCEPTION( set.nb_fixed_objects() == 1, "still one fixed" );
        OPENGEODE_EXCEPTION( set.nb_free_objects() == 2, "two remining free" );
        OPENGEODE_EXCEPTION( set.nb_objects() == 3, "three object at last" );
    }
} // namespace

int main()
{
    try
    {
        geode::Logger::info( "TEST ObjectSet" );
        geode::StochasticLibrary::initialize();

        test_add_and_access();
        test_update_object();
        test_remove_object();
        test_const_access();
        test_fixed_and_free_management();
        test_string_representation();

        geode::Logger::info( "TEST ObjectSet SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}