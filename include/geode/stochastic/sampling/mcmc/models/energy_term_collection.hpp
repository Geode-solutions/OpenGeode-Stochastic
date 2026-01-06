#pragma once

#include <absl/container/flat_hash_map.h>
#include <geode/stochastic/sampling/mcmc/models/components/energy_term.hpp>

namespace geode
{
    template < typename ObjectType >
    class EnergyTermCollection
    {
        OPENGEODE_DISABLE_COPY( EnergyTermCollection );

    public:
        EnergyTermCollection() = default;
        EnergyTermCollection( EnergyTermCollection&& ) noexcept = default;
        ~EnergyTermCollection() = default;

        uuid add_energy_term( std::shared_ptr< EnergyTerm< ObjectType > > term )
        {
            const uuid id = term->id();
            energy_terms_.emplace( id, term );
            for( const uuid& set_id : term->targeted_set_ids() )
            {
                set_to_terms_[set_id].push_back( term );
            }
            return id;
        }

        bool remove_energy_term( const uuid& id )
        {
            auto it = energy_terms_.find( id );
            if( it == energy_terms_.end() )
                return false;

            auto term = it->second;

            for( const uuid& set_id : term->targeted_set_ids() )
            {
                auto vec_it = set_to_terms_.find( set_id );
                if( vec_it == set_to_terms_.end() )
                    continue;

                auto& vec = vec_it->second;
                vec.erase(
                    std::remove( vec.begin(), vec.end(), term ), vec.end() );

                if( vec.empty() )
                    set_to_terms_.erase( vec_it );
            }

            energy_terms_.erase( it );
            return true;
        }

        void clear()
        {
            energy_terms_.clear();
            set_to_terms_.clear();
        }

        [[nodiscard]] index_t size() const
        {
            return energy_terms_.size();
        }

        [[nodiscard]] const EnergyTerm< ObjectType >& get(
            const uuid& id ) const
        {
            auto it = energy_terms_.find( id );
            OPENGEODE_EXCEPTION( it != energy_terms_.end(),
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    id.string() ) );
            return *it->second;
        }

        [[nodiscard]] const absl::flat_hash_map< uuid,
            std::shared_ptr< EnergyTerm< ObjectType > > >&
            all_terms() const
        {
            return energy_terms_;
        }

        [[nodiscard]] const std::vector<
            std::shared_ptr< EnergyTerm< ObjectType > > >&
            terms_for_set( const uuid& set_id ) const
        {
            const auto it = set_to_terms_.find( set_id );
            OPENGEODE_EXCEPTION( it != set_to_terms_.end(),
                "[EnergyTermCollection] - Object Subset (", set_id.string(),
                ") does not have any energy term." );
            return it->second;
        }

        std::string string() const
        {
            auto message = absl::StrCat(
                "EnergyTermCollection: ", energy_terms_.size(), " terms:" );
            for( const auto& [id, term] : energy_terms_ )
            {
                absl::StrAppend( &message, "\n\t --> ", term->string() );
            }
            return message;
        }

    private:
        absl::flat_hash_map< uuid, std::shared_ptr< EnergyTerm< ObjectType > > >
            energy_terms_;

        absl::flat_hash_map< uuid,
            std::vector< std::shared_ptr< EnergyTerm< ObjectType > > > >
            set_to_terms_;
    };

} // namespace geode