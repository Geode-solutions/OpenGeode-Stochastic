// A pattern is collection of spatialized objects:
//  * Basic Object
//  * parametric Object-->marked points
//  * meshes
//  a pattern is the superposition of several subcollection of object of
//  différent type. ( fracture sets)

// a pattern should be easy to modify...
// each object can be added, removed or modified several times.
// we also want to have access to neighbors( fist, second or third circle of
// neighbors) rtree?
//
// a collections of helpers need to be done to evaluates "interactions" between
// objects.
//

#pragma once

#include <vector>

#include <geode/stochastic/configuration/marked_object.hpp>

namespace geode
{
    template < typename Geometry >
    class Configuration
    {
    public:
        void add_object( MarkedObject< Geometry >&& object );
        void change_object( index_t idx, MarkedObject< Geometry >&& object );
        void remove_object( index_t idx );
        index_t size() const;

        const MarkedObject< Geometry >& object( index_t idx ) const;

        std::vector< index_t > object_ids_with_mark( const Mark& mark );

    private:
        std::vector< MarkedObject< Geometry > > objects_;
    };
} // namespace geode

// std::string sample_three( absl::BitGen& gen, double p1, double p2, double p3
// )
//{
//     if( absl::Bernoulli( gen, p1 ) )
//     {
//         return "A";
//     }
//     if( absl::Bernoulli( gen, p2 / ( p2 + p3 ) ) )
//     { // conditional probability
//         return "B";
//     }
//     return "C";
// }