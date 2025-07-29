
#include <geode/stochastic/configuration/configuration.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < typename Geometry >
    void Configuration< Geometry >::add_object( MarkedObject< Geometry > obj )
    {
        objects_.emplace_back( std::move( obj ) );
    }

    template < typename Geometry >
    void Configuration< Geometry >::remove_object( std::size_t idx )
    {
        if( idx < objects_.size() )
        {
            objects_.erase( objects_.begin() + idx );
        }
    }

    template < typename Geometry >
    std::size_t Configuration< Geometry >::size() const
    {
        return objects_.size();
    }

    template < typename Geometry >
    MarkedObject< Geometry >& Configuration< Geometry >::object(
        std::size_t idx )
    {
        return objects_.at( idx );
    }

    template < typename Geometry >
    const MarkedObject< Geometry >& Configuration< Geometry >::object(
        std::size_t idx ) const
    {
        return objects_.at( idx );
    }

    template < typename Geometry >
    std::vector< MarkedObject< Geometry >* >
        Configuration< Geometry >::objects_with_mark( int label )
    {
        std::vector< MarkedObject< Geometry >* > result;
        for( auto& obj : objects_ )
        {
            if( obj.has_mark() && obj.mark()->label == label )
            {
                result.push_back( &obj );
            }
        }
        return result;
    }

    template < typename Geometry >
    std::vector< MarkedObject< Geometry >* >
        Configuration< Geometry >::objects_with_mark()
    {
        std::vector< MarkedObject< Geometry >* > result;
        for( auto& obj : objects_ )
        {
            if( obj.has_mark() )
            {
                result.push_back( &obj );
            }
        }
        return result;
    }

    // Explicit instantiations
    template class opengeode_stochastic_stochastic_api Configuration< Point2D >;
    template class opengeode_stochastic_stochastic_api Configuration< Point3D >;
    template class opengeode_stochastic_stochastic_api
        Configuration< OwnerSegment2D >;
} // namespace geode