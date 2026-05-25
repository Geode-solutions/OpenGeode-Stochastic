/*
 * Copyright (c) 2019 - 2026 Geode-solutions
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

namespace geode::detail
{
    struct EnergyScale
    {
        explicit EnergyScale( double param )
        {
            OpenGeodeStochasticStochasticException::check_exception(
                param >= 0., nullptr, OpenGeodeException::TYPE::data,
                "[EnergyTerm] Model parameter cannot be "
                "negative." );

            if( param >= geode::GLOBAL_EPSILON )
            {
                value_ = -std::log( param ); // store log-space parameter
            }
            // else value = std::nullopt → special case: param == 0
        }

        /// Compute energy contribution for a given statistic multiplier
        [[nodiscard]] double contribution( double multiplier ) const
        {
            if( value_ )
            {
                return value_.value() * multiplier;
            }
            // Hard constraint: param == 0
            return ( multiplier > 0 )
                       ? std::numeric_limits< double >::infinity()
                       : 0.0;
        }

        /// Return original parameter (gamma)
        [[nodiscard]] double parameter() const
        {
            if( value_ )
            {
                return std::exp( -value_.value() );
            }
            return 0.;
        }

    private:
        std::optional< double > value_; // empty if param == 0 (hardcore)
    };
} // namespace geode::detail

namespace geode
{
    template < typename ObjectType >
    class EnergyTerm : public Identifier
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( EnergyTerm );

    public:
        explicit EnergyTerm( std::string_view name,
            double param,
            std::vector< uuid >&& impacted_set_ids,
            const SpatialDomain< ObjectType::dim >& domain )
            : energy_scale_{ param },
              impacted_set_ids_{ std::move( impacted_set_ids ) },
              domain_( domain )
        {
            absl::c_sort( impacted_set_ids_ );
            impacted_set_ids_.erase( std::unique( impacted_set_ids_.begin(),
                                         impacted_set_ids_.end() ),
                impacted_set_ids_.end() );
            impacted_set_ids_.shrink_to_fit();
            IdentifierBuilder builder( *this );
            builder.set_name( name );
        }

        virtual ~EnergyTerm() = default;

        [[nodiscard]] double parameter() const
        {
            return energy_scale_.parameter();
        }

        [[nodiscard]] const std::vector< uuid >& impacted_set_ids() const
        {
            return impacted_set_ids_;
        }

        /// Energy contribution for a given statistic multiplier
        [[nodiscard]] double contribution( double multiplier ) const
        {
            return energy_scale_.contribution( multiplier );
        }

        [[nodiscard]] virtual double total_log(
            const ObjectSets< ObjectType >& state ) const = 0;

        [[nodiscard]] virtual double delta_log_add(
            const ObjectSets< ObjectType >& state,
            const ObjectRef< ObjectType >& new_object ) const = 0;

        [[nodiscard]] virtual double delta_log_remove(
            const ObjectSets< ObjectType >& state,
            const ObjectId& object_id ) const = 0;

        [[nodiscard]] virtual double delta_log_change(
            const ObjectSets< ObjectType >& state,
            const ObjectId& old_object_id,
            const ObjectRef< ObjectType >& new_object ) const = 0;

        [[nodiscard]] virtual double statistic(
            const ObjectSets< ObjectType >& state ) const = 0;

        [[nodiscard]] std::string string() const
        {
            auto message =
                absl::StrCat( "Term : ", name().value_or( id().string() ),
                    "; uuid: ", id().string(),
                    " parameter value: ", energy_scale_.parameter(),
                    " applyied on ", impacted_set_ids_.size(),
                    " object subsets -->" );
            for( const auto& set_id : impacted_set_ids_ )
            {
                absl::StrAppend( &message, "\t", set_id.string() );
            }
            return message;
        }

    protected:
        [[nodiscard]] bool is_impacted_set( const uuid& set_id ) const
        {
            return std::binary_search(
                impacted_set_ids_.begin(), impacted_set_ids_.end(), set_id );
        }

        [[nodiscard]] const SpatialDomain< ObjectType::dim >& domain() const
        {
            return domain_;
        }

        template < typename Func >
        void for_each_object_in_sets( const ObjectSets< ObjectType >& state,
            const std::vector< uuid >& set_ids,
            Func&& do_apply ) const
        {
            for( const auto& set_id : set_ids )
            {
                for( const auto object_ids :
                    state.get_objects_in_set( set_id ) )
                {
                    std::forward< Func >( do_apply )( object_ids );
                }
            }
        }

    private:
        detail::EnergyScale energy_scale_;
        std::vector< uuid > impacted_set_ids_;
        const SpatialDomain< ObjectType::dim >& domain_;
    };
} // namespace geode