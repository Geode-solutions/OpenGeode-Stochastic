/*
 * Copyright (c) 2019 - 2025 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#pragma once

#include <absl/container/flat_hash_set.h>

#include <geode/basic/identifier.hpp>
#include <geode/basic/identifier_builder.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>

#include <geode/stochastic/spatial/spatial_domain.hpp>

#include <optional>

namespace geode
{
    namespace detail
    {
        struct EnergyScale
        {
            explicit EnergyScale( double param )
            {
                OPENGEODE_EXCEPTION( param >= 0.,
                    "[Gibbs energy term] - The model parameter "
                    "cannot be negative." );

                if( param >= geode::GLOBAL_EPSILON )
                {
                    value = -std::log( param ); // store log-space parameter
                }
                // else value = std::nullopt → special case: param == 0
            }

            /// Compute energy contribution for a given statistic multiplier
            double contribution( double multiplier ) const
            {
                if( value )
                {
                    return value.value() * multiplier;
                }
                // Hard constraint: param == 0
                return ( multiplier > 0 )
                           ? std::numeric_limits< double >::infinity()
                           : 0.0;
            }

            /// Return original parameter (gamma)
            double parameter() const
            {
                if( value )
                {
                    return std::exp( -value.value() );
                }
                return 0.;
            }

        private:
            std::optional< double > value; // empty if param == 0 (hardcore)
        };
    } // namespace detail
} // namespace geode

namespace geode
{
    template < typename ObjectType >
    class EnergyTerm : public Identifier
    {
    public:
        explicit EnergyTerm( std::string_view name,
            double param,
            std::vector< uuid > &&targeted_set_ids,
            const SpatialDomain< ObjectType::dim > &domain )
            : energy_scale_{ param },
              targeted_set_ids_{ std::move( targeted_set_ids ) },
              domain_( domain )
        {
            std::sort( targeted_set_ids_.begin(), targeted_set_ids_.end() );
            IdentifierBuilder builder( *this );
            builder.set_name( name );
        }

        virtual ~EnergyTerm() = default;

        double parameter() const
        {
            return energy_scale_.parameter();
        }

        const std::vector< uuid > &targeted_set_ids() const
        {
            return targeted_set_ids_;
        }

        /// Energy contribution for a given statistic multiplier
        double contribution( double multiplier ) const
        {
            return energy_scale_.contribution( multiplier );
        }

        virtual double total_log(
            const ObjectSets< ObjectType > &state ) const = 0;

        virtual double delta_log_add( const ObjectSets< ObjectType > &state,
            const ObjectRef< ObjectType > &new_object ) const = 0;

        virtual double delta_log_remove( const ObjectSets< ObjectType > &state,
            const ObjectId &object_id ) const = 0;

        virtual double delta_log_change( const ObjectSets< ObjectType > &state,
            const ObjectId &old_object_id,
            const ObjectRef< ObjectType > &new_object ) const = 0;

        virtual double statistic(
            const ObjectSets< ObjectType > &state ) const = 0;

        std::string string() const
        {
            auto message =
                absl::StrCat( "Term : ", name(), "; uuid: ", id().string(),
                    " parameter value: ", energy_scale_.parameter(),
                    " applyied on ", targeted_set_ids_.size(),
                    " object subsets -->" );
            for( const auto &set_id : targeted_set_ids_ )
            {
                absl::StrAppend( &message, "\t", set_id.string() );
            }
            return message;
        }

    protected:
        bool is_targeted_set( const uuid &set_id ) const
        {
            return std::binary_search(
                targeted_set_ids_.begin(), targeted_set_ids_.end(), set_id );
        }

        const SpatialDomain< ObjectType::dim > &domain() const
        {
            return domain_;
        }

        //        bool intersects_domain( const ObjectType& obj ) const
        //        {
        //            return SpatialDomainChecker< ObjectType
        //            >::intersects_domain(
        //                domain_, obj );
        //        }

        template < typename Func >
        void for_each_targeted_object(
            const ObjectSets< ObjectType > &state, Func &&do_apply ) const
        {
            for( const auto &targeted_set_id : targeted_set_ids_ )
            {
                for( const auto id :
                    state.get_objects_in_set( targeted_set_id ) )
                {
                    do_apply( id );
                }
            }
        }

    private:
        detail::EnergyScale energy_scale_;
        std::vector< uuid > targeted_set_ids_;
        const SpatialDomain< ObjectType::dim > &domain_;
    };
} // namespace geode