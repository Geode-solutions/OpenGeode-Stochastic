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

        uuid add_energy_term( std::unique_ptr< EnergyTerm< ObjectType > > term )
        {
            const uuid id = term->id();
            EnergyTerm< ObjectType >* ptr = term.get();

            energy_terms_.emplace( id, std::move( term ) );

            for( const uuid& subset_id : ptr->targeted_subset_ids() )
            {
                subset_to_terms_[subset_id].push_back( ptr );
            }
            return id;
        }

        bool remove_energy_term( const uuid& id )
        {
            auto it = energy_terms_.find( id );
            if( it == energy_terms_.end() )
                return false;

            EnergyTerm< ObjectType >* ptr = it->second.get();
            for( const uuid& subset_id : ptr->targeted_subset_ids() )
            {
                auto& vec = subset_to_terms_[subset_id];
                vec.erase(
                    std::remove( vec.begin(), vec.end(), ptr ), vec.end() );
                if( vec.empty() )
                    subset_to_terms_.erase( subset_id );
            }

            energy_terms_.erase( it );
            return true;
        }

        void clear()
        {
            energy_terms_.clear();
            subset_to_terms_.clear();
        }

        [[nodiscard]] index_t size() const
        {
            return energy_terms_.size();
        }

        [[nodiscard]] EnergyTerm< ObjectType >* get( const uuid& id ) const
        {
            auto it = energy_terms_.find( id );
            OPENGEODE_EXCEPTION( it != energy_terms_.end(),
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    id.string() ) );
            return it->second.get();
        }

        [[nodiscard]] std::vector< EnergyTerm< ObjectType >* > terms_for_subset(
            const uuid& subset_id ) const
        {
            auto it = subset_to_terms_.find( subset_id );
            return it != subset_to_terms_.end()
                       ? it->second
                       : std::vector< EnergyTerm< ObjectType >* >{};
        }

        [[nodiscard]] std::vector< EnergyTerm< ObjectType >* > all_terms() const
        {
            std::vector< EnergyTerm< ObjectType >* > result;
            result.reserve( energy_terms_.size() );
            for( const auto& [id, term] : energy_terms_ )
            {
                result.push_back( term.get() );
            }
            return result;
        }

        std::string string() const
        {
            auto message = absl::StrCat(
                "Gibbs Energy: ", energy_terms_.size(), " terms:" );
            // we should find a way to iterate in a ordered way.
            for( const auto energy_term : all_terms() )
            {
                absl::StrAppend( &message, "\n\t --> ", energy_term->string() );
            }
            return message;
        }

    private:
        absl::flat_hash_map< uuid, std::unique_ptr< EnergyTerm< ObjectType > > >
            energy_terms_;
        absl::flat_hash_map< uuid, std::vector< EnergyTerm< ObjectType >* > >
            subset_to_terms_;
    };
} // namespace geode