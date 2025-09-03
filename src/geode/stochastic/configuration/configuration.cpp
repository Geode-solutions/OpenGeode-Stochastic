
#include <geode/stochastic/configuration/configuration.hpp>
#include <geode/stochastic/configuration/object_helpers.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < typename Object >
    const Object& Configuration< Object >::get_object( ObjectId id ) const
    {
        const auto& group = get_group( id.group );
        OPENGEODE_EXCEPTION( id.object < group.size(),
            "[Configuration]- object index out of range." );
        return group[id.object];
    }

    template < typename Object >
    const std::vector< Object >& Configuration< Object >::get_group(
        const GroupId& group_id ) const
    {
        auto it = groups_.find( group_id );
        OPENGEODE_EXCEPTION( it != groups_.end(), "[Configuration] - group (",
            group_id.value, ") is not defined." );
        return it->second;
    }

    template < typename Object >
    std::vector< ObjectId > Configuration< Object >::get_all_object() const
    {
        std::vector< ObjectId > result;
        result.reserve( nb_objects() );
        for( const auto& [group_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                result.push_back( { obj_id, group_id } );
            }
        }
        return result;
    }

    template < typename Object >
    index_t Configuration< Object >::nb_groups() const
    {
        return groups_.size();
    }

    template < typename Object >
    index_t Configuration< Object >::nb_objects_in_group(
        const GroupId& group_id ) const
    {
        return get_group( group_id ).size();
    }

    template < typename Object >
    index_t Configuration< Object >::nb_objects() const
    {
        index_t nb_objects{ 0 };
        for( const auto& [group_id, objs] : groups_ )
        {
            geode_unused( group_id );
            nb_objects += objs.size();
        }
        return nb_objects;
    }

    template < typename Object >
    ObjectId Configuration< Object >::add_object(
        Object&& object, GroupId group_id )
    {
        auto& group = groups_[group_id]; // creates group if missing
        ObjectId new_object_id{ static_cast< index_t >( group.size() ),
            group_id };
        group.push_back( std::move( object ) );
        //  tree_.add( object_bounding_box( object ), object_id );
        return new_object_id;
    }

    template < typename Object >
    void Configuration< Object >::update_object( ObjectId id, Object&& object )
    {
        auto& group = get_group( id.group );
        OPENGEODE_EXCEPTION( id.object < group.size(),
            "[Configuration]- index of object to update out of range." );
        group[id.object] = std::move( object );
    }

    template < typename Object >
    void Configuration< Object >::remove_object( ObjectId id )
    {
        auto& group = get_group( id.group );
        OPENGEODE_EXCEPTION( id.object < group.size(),
            "[Configuration]- index of object to remove out of range." );
        group[id.object] = std::move( group.back() );
        group.pop_back();
    }
    template < typename Object >
    std::vector< ObjectId > Configuration< Object >::neighbors(
        const ObjectId object_id, double searching_distance ) const
    {
        geode_unused( searching_distance );
        std::vector< ObjectId > result;
        for( const auto& [group_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                const ObjectId cur_object_id{ obj_id, group_id };
                if( cur_object_id == object_id )
                {
                    continue;
                }
                result.push_back( cur_object_id );
            }
        }
        return result;
    }

    template < typename Object >
    std::vector< ObjectId > Configuration< Object >::neighbors(
        const Object& object, double searching_distance ) const
    {
        geode_unused( searching_distance );
        geode_unused( object );
        std::vector< ObjectId > result;
        for( const auto& [group_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                const ObjectId cur_object_id{ obj_id, group_id };
                result.push_back( cur_object_id );
            }
        }
        return result;
    }

    template < typename Object >
    std::vector< Object >& Configuration< Object >::get_group(
        const GroupId& group_id )
    {
        return const_cast< std::vector< Object >& >(
            static_cast< const Configuration* >( this )->get_group(
                group_id ) );
    }

    template class opengeode_stochastic_stochastic_api Configuration< Point2D >;
    template class opengeode_stochastic_stochastic_api Configuration< Point3D >;
    template class opengeode_stochastic_stochastic_api
        Configuration< OwnerSegment2D >;
} // namespace geode