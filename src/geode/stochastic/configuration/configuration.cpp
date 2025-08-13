
#include <geode/stochastic/configuration/configuration.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < typename Geometry >
    void Configuration< Geometry >::add_object(
        MarkedObject< Geometry >&& object )
    {
        objects_.emplace_back( std::move( object ) );
    }

    template < typename Geometry >
    void Configuration< Geometry >::change_object(
        index_t idx, MarkedObject< Geometry >&& object )
    {
        OPENGEODE_EXCEPTION( idx < objects_.size(),
            "[Configuration]- index to change out of range;" );
        objects_[idx] = std::move( object );
    }

    template < typename Geometry >
    void Configuration< Geometry >::remove_object( index_t idx )
    {
        OPENGEODE_EXCEPTION( idx < objects_.size(),
            "[Configuration]- index to remove out of range;" );
        objects_[idx] = std::move( objects_.back() );
        objects_.pop_back();
    }

    template < typename Geometry >
    index_t Configuration< Geometry >::size() const
    {
        return objects_.size();
    }

    template < typename Geometry >
    const MarkedObject< Geometry >& Configuration< Geometry >::object(
        index_t idx ) const
    {
        return objects_.at( idx );
    }

    template < typename Geometry >
    const MarkedObject< Geometry >& Configuration< Geometry >::operator[](
        index_t idx ) const
    {
        return objects_[idx];
    }

    template < typename Geometry >
    std::vector< index_t > Configuration< Geometry >::object_ids_with_mark(
        const Mark& mark )
    {
        std::vector< index_t > result;
        result.reserve( objects_.size() );
        for( const auto idx : Range{ objects_.size() } )
        {
            const auto& obj = objects_[idx];
            if( obj.has_mark() && obj.mark() == mark )
            {
                result.emplace_back( idx );
            }
        }
        return result;
    }

    template class opengeode_stochastic_stochastic_api Configuration< Point2D >;
    template class opengeode_stochastic_stochastic_api Configuration< Point3D >;
    template class opengeode_stochastic_stochastic_api
        Configuration< OwnerSegment2D >;
} // namespace geode