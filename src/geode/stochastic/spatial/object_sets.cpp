
#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < typename Type >
    const ObjectSet< Type >& ObjectSets< Type >::get_set(
        const uuid& set_id ) const
    {
        auto it = uuid_to_index_.find( set_id );
        OPENGEODE_EXCEPTION( it != uuid_to_index_.end(),
            "[ObjectSet] - group (", set_id.string(), ") is not defined." );
        return sets_[it->second];
    }

    template < typename Type >
    const Type& ObjectSets< Type >::get_object( const ObjectId& object ) const
    {
        auto& set = get_set( object.set_id );
        if( object.fixed )
        {
            return set.get_fixed_object( object.index );
        }
        return set.get_free_object( object.index );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::get_all_object() const
    {
        std::vector< ObjectId > result;
        result.reserve( nb_objects() );
        for( const auto& objs : sets_ )
        {
            auto set_id = objs.id();
            for( const auto obj_id : geode::Range{ objs.nb_fixed_objects() } )
            {
                result.push_back( { obj_id, true, set_id } );
            }
            for( const auto obj_id : geode::Range{ objs.nb_free_objects() } )
            {
                result.push_back( { obj_id, false, set_id } );
            }
        }
        return result;
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::get_objects_in_set(
        const uuid& set_id ) const
    {
        const auto& set = get_set( set_id );
        std::vector< ObjectId > result;
        result.reserve( set.nb_objects() );
        for( const auto obj_id : geode::Range{ set.nb_fixed_objects() } )
        {
            result.push_back( { obj_id, true, set_id } );
        }
        for( const auto obj_id : geode::Range{ set.nb_free_objects() } )
        {
            result.push_back( { obj_id, false, set_id } );
        }
        return result;
    }

    template < typename Type >
    index_t ObjectSets< Type >::nb_sets() const
    {
        return sets_.size();
    }

    template < typename Type >
    index_t ObjectSets< Type >::nb_objects_in_set( const uuid& set_id ) const
    {
        return get_set( set_id ).nb_objects();
    }

    template < typename Type >
    index_t ObjectSets< Type >::nb_objects() const
    {
        index_t nb_objects{ 0 };
        for( const auto& objs : sets_ )
        {
            nb_objects += objs.nb_objects();
        }
        return nb_objects;
    }

    template < typename Type >
    uuid ObjectSets< Type >::add_set( std::string_view name )
    {
        auto set_index = sets_.size();
        auto& new_set = sets_.emplace_back( ObjectSet< Type >{} );
        const auto set_uuid = new_set.id();

        new_set.set_name( name );
        auto [it_set_name, set_uuid_inserted] =
            name_to_uuid_.emplace( name, set_uuid );
        OPENGEODE_EXCEPTION(
            set_uuid_inserted, absl::StrCat( "[ObjectSet]- group named ", name,
                                   " already exists." ) );

        auto [it_set_uuid, set_index_inserted] =
            uuid_to_index_.emplace( set_uuid, set_index );
        OPENGEODE_EXCEPTION( set_index_inserted, "[ObjectSet]- group (",
            set_uuid.string(), ") already exists." );
        return set_uuid;
    }

    template < typename Type >
    ObjectId ObjectSets< Type >::add_object(
        Type&& object, const uuid& set_id, bool fixed )
    {
        auto& set = get_set( set_id );
        auto object_box = object_bounding_box( object );

        if( fixed )
        {
            auto new_fixed_id = set.add_fixed_object( std::move( object ) );
            ObjectId fixed_object_id{ new_fixed_id, true, set_id };
            neighborhood_.add( object_box, fixed_object_id );
            return fixed_object_id;
        }
        auto new_free_id = set.add_free_object( std::move( object ) );
        ObjectId free_object_id{ new_free_id, false, set_id };
        neighborhood_.add( object_box, free_object_id );
        return free_object_id;
    }

    template < typename Type >
    void ObjectSets< Type >::update_free_object(
        const ObjectId& old_object_id, Type&& new_object )
    {
        OPENGEODE_EXCEPTION(
            !old_object_id.fixed, "[ObjectSet]- cannot modify fixed object." );
        auto& set = get_set( old_object_id.set_id );
        OPENGEODE_EXCEPTION( old_object_id.index < set.nb_objects(),
            "[ObjectSet]- index of object to update out of range." );
        auto old_box = object_bounding_box( get_object( old_object_id ) );
        auto new_box = object_bounding_box( new_object );
        neighborhood_.update( old_box, new_box, old_object_id );
        set.update_free_object( old_object_id.index, std::move( new_object ) );
    }

    template < typename Type >
    void ObjectSets< Type >::remove_free_object( const ObjectId& object_id )
    {
        auto& set = get_set( object_id.set_id );
        OPENGEODE_EXCEPTION(
            !object_id.fixed, "[ObjectSet]- Cannot remove fixed object." );
        const auto& obj_to_remove = get_object( object_id );
        neighborhood_.remove( object_bounding_box( obj_to_remove ), object_id );

        ObjectId last_free_id{ static_cast< geode::index_t >(
                                   set.nb_free_objects() - 1 ),
            false, object_id.set_id };
        if( object_id != last_free_id )
        {
            const auto& last_free_obj =
                set.get_free_object( last_free_id.index );
            auto box_to_move = object_bounding_box( last_free_obj );
            neighborhood_.update( box_to_move, last_free_id, object_id );
        }
        set.remove_free_object( object_id.index );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::neighbors( const Type& object,
        const std::vector< uuid >& targeted_set_ids,
        double searching_distance,
        std::optional< ObjectId > excluded_id ) const
    {
        auto box = object_bounding_box( object );
        box.extends( searching_distance * 2. );
        return neighborhood_.get_all_neighbor_ids(
            box, targeted_set_ids, excluded_id );
    }

    template < typename Type >
    ObjectSet< Type >& ObjectSets< Type >::get_set( const uuid& set_id )
    {
        return const_cast< ObjectSet< Type >& >(
            static_cast< const ObjectSets* >( this )->get_set( set_id ) );
    }

    template < typename Type >
    uuid ObjectSets< Type >::get_set_uuid( const std::string_view name ) const
    {
        if( auto set_uuid = name_to_uuid_.find( name );
            set_uuid != name_to_uuid_.end() )
        {
            return set_uuid->second;
        }
        throw OpenGeodeException(
            "[ObjectSets] ObjectSet uuid accessor - group named ", name,
            " does not exist." );
    }

    template < typename Type >
    std::vector< uuid > ObjectSets< Type >::get_set_uuids(
        const std::vector< std::string >& set_names ) const
    {
        std::vector< geode::uuid > uuids;
        uuids.reserve( set_names.size() );

        for( const auto& name : set_names )
        {
            uuids.emplace_back( get_set_uuid( name ) );
        }
        return uuids;
    }

    template < typename Type >
    std::vector< std::pair< uuid, uuid > >
        ObjectSets< Type >::get_set_uuid_pairs(
            const std::vector< std::pair< std::string, std::string > >&
                set_name_pairs ) const
    {
        std::vector< std::pair< uuid, uuid > > result;
        result.reserve( set_name_pairs.size() );

        for( const auto& [name1, name2] : set_name_pairs )
        {
            result.emplace_back( get_set_uuid( name1 ), get_set_uuid( name2 ) );
        }
        return result;
    }

    template < typename Type >
    std::string ObjectSets< Type >::string() const
    {
        auto message = absl::StrCat( "ObjectSets with ", nb_objects(),
            " objects in, ", nb_sets(), " sets" );
        for( const auto& objs : sets_ )
        {
            absl::StrAppend( &message, "\n\t --> ", objs.string() );
        }
        return message;
    }

    template class opengeode_stochastic_stochastic_api ObjectSets< Point2D >;
    template class opengeode_stochastic_stochastic_api ObjectSets< Point3D >;
    template class opengeode_stochastic_stochastic_api
        ObjectSets< OwnerSegment2D >;
} // namespace geode