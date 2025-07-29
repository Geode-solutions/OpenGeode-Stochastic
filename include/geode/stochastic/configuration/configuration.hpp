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
#include <geode/stochastic/configuration/marked_object.hpp>
#include <vector>
namespace geode
{
    template < typename Geometry >
    class Configuration
    {
    public:
        void add_object( MarkedObject< Geometry > obj );
        void remove_object( std::size_t idx );
        std::size_t size() const;

        MarkedObject< Geometry >& object( std::size_t idx );
        const MarkedObject< Geometry >& object( std::size_t idx ) const;

        std::vector< MarkedObject< Geometry >* > objects_with_mark( int label );
        std::vector< MarkedObject< Geometry >* > objects_with_mark();

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