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
#include <geode/stochastic/common.hpp>
#include <geode/stochastic/spatial/object_set.hpp>
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
    //    struct EnergyTermDescription
    //    {
    //        geode::uuid id;
    //        std::string name;
    //        std::string type;
    //        double parameter_value;
    //        std::optional< uuid > targeted_subset_id{};
    //    }
    //
    //    struct StatisticalDescription
    //    {
    //        std::string label;
    //        double value;
    //        std::optional< uuid > targeted_subset_id{};
    //    };

    template < typename ObjectType >
    class EnergyTerm
    {
    public:
        explicit EnergyTerm( std::string_view name,
            double param,
            absl::flat_hash_set< uuid > targeted_subset_ids )
            : name_{ name },
              energy_scale_{ param },
              targeted_subset_ids_{ std::move( targeted_subset_ids ) }
        {
        }

        virtual ~EnergyTerm() = default;

        const uuid& id() const
        {
            return energy_term_id_;
        }

        std::string_view name() const
        {
            return name_;
        }

        double parameter() const
        {
            return energy_scale_.parameter();
        }

        const absl::flat_hash_set< uuid >& targeted_subset_ids() const
        {
            return targeted_subset_ids_;
        }

        /// Energy contribution for a given statistic multiplier
        double contribution( double multiplier ) const
        {
            return energy_scale_.contribution( multiplier );
        }

        virtual double total_log(
            const ObjectSet< ObjectType >& state ) const = 0;

        virtual double delta_log_add( const ObjectSet< ObjectType >& state,
            const ObjectRef< ObjectType >& new_object ) const = 0;

        virtual double delta_log_remove( const ObjectSet< ObjectType >& state,
            const ObjectId& object_id ) const = 0;

        virtual double delta_log_change( const ObjectSet< ObjectType >& state,
            const ObjectId& old_object_id,
            const ObjectRef< ObjectType >& new_object ) const = 0;

        virtual double statistic(
            const ObjectSet< ObjectType >& state ) const = 0;

        std::string string() const
        {
            auto message =
                absl::StrCat( "Term : ", name(), "; uuid: ", id().string(),
                    " parameter value: ", energy_scale_.parameter(),
                    " applyied on ", targeted_subset_ids_.size(),
                    " object subsets -->" );
            for( const auto& subset_uuid : targeted_subset_ids_ )
            {
                absl::StrAppend( &message, "\t", subset_uuid.string() );
            }
            return message;
        }

    protected:
        bool is_targeted_subset( const uuid& subset_id ) const
        {
            return targeted_subset_ids_.find( subset_id )
                   != targeted_subset_ids_.end();
        }

        template < typename Func >
        void for_each_targeted_object(
            const ObjectSet< ObjectType >& state, Func&& do_apply ) const
        {
            for( const auto& targeted_subset_id : targeted_subset_ids_ )
            {
                for( const auto id : geode::Range{
                         state.nb_objects_in_subset( targeted_subset_id ) } )
                {
                    do_apply( ObjectId{ id, targeted_subset_id } );
                }
            }
        }

    private:
        std::string name_;
        detail::EnergyScale energy_scale_;

        absl::flat_hash_set< uuid > targeted_subset_ids_;
        uuid energy_term_id_{};
    };
} // namespace geode