
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
        auto it = sets_.find( set_id );
        OPENGEODE_EXCEPTION( it != sets_.end(), "[ObjectSet] - group (",
            set_id.string(), ") is not defined." );
        return it->second;
    }

    template < typename Type >
    const Type& ObjectSets< Type >::get_object( const ObjectId& object ) const
    {
        auto& set = get_set( object.set_id );
        OPENGEODE_EXCEPTION( object.index < set.nb_objects(),
            "[ObjectSet]- object index out of range." );
        return set.get_object( object.index );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::get_all_object() const
    {
        std::vector< ObjectId > result;
        result.reserve( nb_objects() );
        for( const auto& [set_id, objs] : sets_ )
        {
            for( const auto obj_id : geode::Range{ objs.nb_objects() } )
            {
                result.push_back( { obj_id, set_id } );
            }
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
        for( const auto& [set_id, objs] : sets_ )
        {
            geode_unused( set_id );
            nb_objects += objs.nb_objects();
        }
        return nb_objects;
    }

    template < typename Type >
    uuid ObjectSets< Type >::add_set( std::string_view name )
    {
        ObjectSet< Type > new_set;
        new_set.set_name( name );
        const auto new_set_id = new_set.id();
        auto [it, inserted] = sets_.emplace( new_set_id, std::move( new_set ) );
        OPENGEODE_EXCEPTION( inserted, "[ObjectSet]- group (",
            new_set_id.string(), ") already exists." );
        return new_set_id;
    }

    template < typename Type >
    ObjectId ObjectSets< Type >::add_fixed_object(
        Type&& object, const uuid& set_id )
    {
        auto& set = get_set( set_id );
        auto new_fixed_id = set.add_fixed_object( std::move( object ) );

        ObjectId fixed_object_id{ new_fixed_id, set_id };
        auto fixed_object = get_object( fixed_object_id );
        auto fixed_object_box = object_bounding_box( fixed_object );

        ObjectId last_object_id{
            static_cast< geode::index_t >( set.nb_objects() - 1 ), set_id
        };
        auto free_object = get_object( last_object_id );
        auto free_object_box = object_bounding_box( free_object );

        neighborhood_.update(
            free_object_box, fixed_object_box, fixed_object_id );
        neighborhood_.add( free_object_box, last_object_id );

        return fixed_object_id;
    }

    template < typename Type >
    ObjectId ObjectSets< Type >::add_free_object(
        Type&& object, const uuid& set_id )
    {
        auto& set = get_set( set_id );
        ObjectId new_object_id{ static_cast< index_t >( set.nb_objects() ),
            set_id };
        neighborhood_.add( object_bounding_box( object ), new_object_id );
        set.add_free_object( std::move( object ) );
        return new_object_id;
    }

    template < typename Type >
    void ObjectSets< Type >::update_free_object(
        const ObjectId& old_object_id, Type&& new_object )
    {
        auto& set = get_set( old_object_id.set_id );
        OPENGEODE_EXCEPTION( old_object_id.index < set.nb_objects(),
            "[ObjectSet]- index of object to update out of range." );
        auto old_box =
            object_bounding_box( set.get_object( old_object_id.index ) );
        auto new_box = object_bounding_box( new_object );
        neighborhood_.update( old_box, new_box, old_object_id );
        set.update_free_object( old_object_id.index, std::move( new_object ) );
    }

    template < typename Type >
    void ObjectSets< Type >::remove_free_object( const ObjectId& object_id )
    {
        auto& set = get_set( object_id.set_id );
        OPENGEODE_EXCEPTION( !set.is_fixed( object_id.index ),
            "[ObjectSet]- Object to remove is fixed." );
        remove_object( object_id );
    }

    template < typename Type >
    void ObjectSets< Type >::remove_object( const ObjectId& object_id )
    {
        auto& set = get_set( object_id.set_id );
        OPENGEODE_EXCEPTION( object_id.index < set.nb_objects(),
            "[ObjectSet]- index of object to remove out of range." );
        const auto& obj_to_remove = set.get_object( object_id.index );
        neighborhood_.remove( object_bounding_box( obj_to_remove ), object_id );

        ObjectId last_id{ static_cast< geode::index_t >( set.nb_objects() - 1 ),
            object_id.set_id };
        if( object_id != last_id )
        {
            const auto& last_obj = set.get_object( last_id.index );
            auto box_to_move = object_bounding_box( last_obj );
            neighborhood_.update( box_to_move, last_id, object_id );
        }
        set.remove_object( object_id.index );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::neighbors(
        const ObjectId& object_id,
        const std::vector< uuid >& targeted_set_ids,
        double searching_distance ) const
    {
        auto box = object_bounding_box( get_object( object_id ) );
        box.extends( searching_distance * 2. );
        return neighborhood_.get_all_neighbor_ids(
            box, targeted_set_ids, object_id );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::neighbors( const Type& object,
        const std::vector< uuid >& targeted_set_ids,
        double searching_distance ) const
    {
        auto box = object_bounding_box( object );
        box.extends( searching_distance * 2. );
        return neighborhood_.get_all_neighbor_ids(
            box, targeted_set_ids, std::nullopt );
    }

    template < typename Type >
    ObjectSet< Type >& ObjectSets< Type >::get_set( const uuid& set_id )
    {
        return const_cast< ObjectSet< Type >& >(
            static_cast< const ObjectSets* >( this )->get_set( set_id ) );
    }

    template < typename Type >
    std::string ObjectSets< Type >::string() const
    {
        auto message = absl::StrCat( "ObjectSets with ", nb_objects(),
            " objects in, ", nb_sets(), " sets" );
        for( const auto& [set_id, objs] : sets_ )
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