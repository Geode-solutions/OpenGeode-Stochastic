
#include <geode/stochastic/configuration/configuration.hpp>
#include <geode/stochastic/configuration/object_helpers.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/point.hpp>

namespace geode
{
    template < typename Type >
    const std::vector< Type >& ObjectSet< Type >::get_subset(
        const uuid& subset_id ) const
    {
        auto it = groups_.find( subset_id );
        OPENGEODE_EXCEPTION( it != groups_.end(), "[ObjectSet] - group (",
            subset_id.string(), ") is not defined." );
        return it->second;
    }

    template < typename Type >
    const Type& ObjectSet< Type >::get_object( const ObjectId& object ) const
    {
        auto& subset = get_subset( object.subset );
        OPENGEODE_EXCEPTION( object.index < subset.size(),
            "[ObjectSet]- object index out of range." );
        return subset[object.index];
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSet< Type >::get_all_object() const
    {
        std::vector< ObjectId > result;
        result.reserve( nb_objects() );
        for( const auto& [subset_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                result.push_back( { obj_id, subset_id } );
            }
        }
        return result;
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_subsets() const
    {
        return groups_.size();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_objects_in_subset(
        const uuid& subset_id ) const
    {
        return get_subset( subset_id ).size();
    }

    template < typename Type >
    index_t ObjectSet< Type >::nb_objects() const
    {
        index_t nb_objects{ 0 };
        for( const auto& [subset_id, objs] : groups_ )
        {
            geode_unused( subset_id );
            nb_objects += objs.size();
        }
        return nb_objects;
    }

    template < typename Type >
    void ObjectSet< Type >::add_subset( const uuid& subset_id )
    {
        auto [it, inserted] =
            groups_.emplace( subset_id, std::vector< Type >{} );
        OPENGEODE_EXCEPTION( inserted, "[ObjectSet]- group (",
            subset_id.string(), ") already exists." );
    }

    template < typename Type >
    ObjectId ObjectSet< Type >::add_object(
        Type&& object, const uuid& subset_id )
    {
        auto& group = get_subset( subset_id );
        ObjectId new_object_id{ static_cast< index_t >( group.size() ),
            subset_id };
        group.push_back( std::move( object ) );
        //  tree_.add( object_bounding_box( object ), object_id );
        return new_object_id;
    }

    template < typename Type >
    void ObjectSet< Type >::update_object(
        const ObjectId& old_object, Type&& new_object )
    {
        auto& subset = get_subset( old_object.subset );
        OPENGEODE_EXCEPTION( old_object.index < subset.size(),
            "[ObjectSet]- index of object to update out of range." );
        subset[old_object.index] = std::move( new_object );
    }

    template < typename Type >
    void ObjectSet< Type >::remove_object( const ObjectId& object )
    {
        auto& subset = get_subset( object.subset );
        OPENGEODE_EXCEPTION( object.index < subset.size(),
            "[ObjectSet]- index of object to remove out of range." );
        subset[object.index] = std::move( subset.back() );
        subset.pop_back();
    }
    template < typename Type >
    std::vector< ObjectId > ObjectSet< Type >::neighbors(
        const ObjectId& object_id, double searching_distance ) const
    {
        geode_unused( searching_distance );
        std::vector< ObjectId > result;
        for( const auto& [subset_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                const ObjectId cur_object_id{ obj_id, subset_id };
                if( cur_object_id == object_id )
                {
                    continue;
                }
                result.push_back( cur_object_id );
            }
        }
        return result;
    }

    template < typename Type >
    std::vector< ObjectId > ObjectSet< Type >::neighbors(
        const Type& object, double searching_distance ) const
    {
        geode_unused( searching_distance );
        geode_unused( object );
        std::vector< ObjectId > result;
        for( const auto& [subset_id, objs] : groups_ )
        {
            for( const auto obj_id : geode::Range{ objs.size() } )
            {
                const ObjectId cur_object_id{ obj_id, subset_id };
                result.push_back( cur_object_id );
            }
        }
        return result;
    }

    template < typename Type >
    std::vector< Type >& ObjectSet< Type >::get_subset( const uuid& subset_id )
    {
        return const_cast< std::vector< Type >& >(
            static_cast< const ObjectSet* >( this )->get_subset( subset_id ) );
    }

    template class opengeode_stochastic_stochastic_api ObjectSet< Point2D >;
    template class opengeode_stochastic_stochastic_api ObjectSet< Point3D >;
    template class opengeode_stochastic_stochastic_api
        ObjectSet< OwnerSegment2D >;
} // namespace geode