
#include <geode/stochastic/spatial/object_helpers.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>
namespace
{
    template < geode::index_t dimension >
    geode::Vector< dimension > enlarge_vector( double distance );

    template <>
    geode::Vector2D enlarge_vector< 2 >( double distance )
    {
        return geode::Vector2D{ { distance, distance } };
    }

    template <>
    geode::Vector3D enlarge_vector< 3 >( double distance )
    {
        return geode::Vector3D{ { distance, distance, distance } };
    }
    template < geode::index_t dimension >
    void enlarge_box(
        geode::BoundingBox< dimension >& box, double search_distance )
    {
        auto vec = enlarge_vector< dimension >( search_distance );
        box.add_point( box.min() - vec );
        box.add_point( box.max() + vec );
    }
} // namespace
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
        OPENGEODE_EXCEPTION( object.index < set.size(),
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
            for( const auto obj_id : geode::Range{ objs.size() } )
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
        return get_set( set_id ).size();
    }

    template < typename Type >
    index_t ObjectSets< Type >::nb_objects() const
    {
        index_t nb_objects{ 0 };
        for( const auto& [set_id, objs] : sets_ )
        {
            geode_unused( set_id );
            nb_objects += objs.size();
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
    ObjectId ObjectSets< Type >::add_object( Type&& object, const uuid& set_id )
    {
        auto& set = get_set( set_id );
        ObjectId new_object_id{ static_cast< index_t >( set.size() ), set_id };
        neighborhood_.add( object_bounding_box( object ), new_object_id );
        set.add_object( std::move( object ) );
        return new_object_id;
    }

    template < typename Type >
    void ObjectSets< Type >::update_object(
        const ObjectId& old_object_id, Type&& new_object )
    {
        auto& set = get_set( old_object_id.set_id );
        OPENGEODE_EXCEPTION( old_object_id.index < set.size(),
            "[ObjectSet]- index of object to update out of range." );
        auto old_box =
            object_bounding_box( set.get_object( old_object_id.index ) );
        auto new_box = object_bounding_box( new_object );
        neighborhood_.update( old_box, new_box, old_object_id );
        set.update_object( old_object_id.index, std::move( new_object ) );
    }

    template < typename Type >
    void ObjectSets< Type >::remove_object( const ObjectId& object_id )
    {
        auto& set = get_set( object_id.set_id );
        OPENGEODE_EXCEPTION( object_id.index < set.size(),
            "[ObjectSet]- index of object to remove out of range." );
        const auto& obj_to_remove = set.get_object( object_id.index );
        neighborhood_.remove( object_bounding_box( obj_to_remove ), object_id );

        ObjectId last_id{ static_cast< geode::index_t >( set.size() - 1 ),
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
        const ObjectId& object_id, double searching_distance ) const
    {
        auto box = object_bounding_box( get_object( object_id ) );
        enlarge_box< Type::dim >( box, searching_distance );
        return neighborhood_.get_all_neighbor_ids( box, object_id );
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSets< Type >::neighbors(
        const Type& object, double searching_distance ) const
    {
        auto box = object_bounding_box( object );
        enlarge_box< Type::dim >( box, searching_distance );
        return neighborhood_.get_all_neighbor_ids( box, std::nullopt );
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
        auto message =
            absl::StrCat( "ObjectSets with ", nb_sets(), " ObjectSet:" );
        for( const auto& [set_id, objs] : sets_ )
        {
            absl::StrAppend( &message, "\n\t - set uuid: ", set_id.string(),
                "--> ", objs.string() );
        }
        return message;
    }

    template class opengeode_stochastic_stochastic_api ObjectSets< Point2D >;
    template class opengeode_stochastic_stochastic_api ObjectSets< Point3D >;
    template class opengeode_stochastic_stochastic_api
        ObjectSets< OwnerSegment2D >;
} // namespace geode