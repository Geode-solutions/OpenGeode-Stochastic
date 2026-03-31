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

#include <geode/basic/common.hpp>
#include <geode/stochastic/models/model_configuration.hpp>

namespace
{
    struct SetDescription
    {
        std::string name;
        geode::index_t number_of_object;
        std::vector< geode::Point2D > fixed_object;
    };
} // namespace

int main()
{
    try
    {
        geode::StochasticLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );

        std::string set_name1 = "set01";
        std::string set_name2 = "set02";
        std::string set_name3 = "set03";
        std::string set_name4 = "set04";

        geode::ObjectSets< geode::Point2D > object_sets;
        const auto set_id1 = object_sets.add_set( set_name1 );
        const auto set_id2 = object_sets.add_set( set_name2 );
        const auto set_id3 = object_sets.add_set( set_name3 );
        const auto set_id4 = object_sets.add_set( set_name4 );

        geode::ModelConfig config;
        config.terms.push_back(
            geode::DensityTermConfig{ "density_set1", { set_name1 }, 1.0 } );
        config.terms.push_back(
            geode::DensityTermConfig{ "density_set2", { set_name2 }, 1.0 } );
        config.terms.push_back(
            geode::DensityTermConfig{ "density_set3", { set_name3 }, 1.0 } );
        config.terms.push_back(
            geode::DensityTermConfig{ "density_set4", { set_name4 }, 1.0 } );
        config.terms.push_back( geode::PairwiseTermConfig{
            "eclidiant_set2_3_4", { set_name4, set_name2, set_name3 }, 0.5 } );

        geode::BoundingBox2D box;
        box.add_point( geode::Point2D{ { 0.0, 0.0 } } );
        box.add_point( geode::Point2D{ { 10.0, 10.0 } } );
        geode::SpatialDomain domain( box, 0. );

        auto collection = geode::build_energy_term_collection< geode::Point2D >(
            config, object_sets, domain );
        geode::Logger::info( collection.size() );
        OPENGEODE_EXCEPTION( collection.size() == 5, "Collection not created" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}