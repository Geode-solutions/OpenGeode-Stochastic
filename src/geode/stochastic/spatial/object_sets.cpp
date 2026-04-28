
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
        OpenGeodeStochasticStochasticException::check( it != sets_.end(),
            nullptr, OpenGeodeException::TYPE::data, "[ObjectSet] - group (",
            set_id.string(), ") is not defined." );
        return it->second;
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
        for( const auto& [set_id, objs] : sets_ )
        {
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
        OpenGeodeStochasticStochasticException::check( inserted, new_set_id,
            OpenGeodeException::TYPE::data, "[ObjectSet]- group (",
            new_set_id.string(), ") already exists." );
        return new_set_id;
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
        OpenGeodeStochasticStochasticException::check( !old_object_id.fixed,
            nullptr, OpenGeodeException::TYPE::data,
            "[ObjectSet]- cannot modify fixed object." );
        auto& set = get_set( old_object_id.set_id );
        OpenGeodeStochasticStochasticException::check(
            old_object_id.index < set.nb_objects(), nullptr,
            OpenGeodeException::TYPE::data,
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
        OpenGeodeStochasticStochasticException::check( !object_id.fixed,
            nullptr, OpenGeodeException::TYPE::data,
            "[ObjectSet]- Cannot remove fixed object." );
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