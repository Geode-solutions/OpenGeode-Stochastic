#pragma once

#include <absl/container/flat_hash_map.h>
#include <geode/stochastic/models/energy_terms/energy_term.hpp>

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

        EnergyTermCollection& operator=( EnergyTermCollection&& ) = default;

        [[nodiscard]] uuid add_energy_term(
            std::unique_ptr< EnergyTerm< ObjectType > >&& term )
        {
            const uuid term_uuid = term->id();
            auto term_idx = energy_terms_.size();
            energy_terms_.emplace_back( std::move( term ) );
            energy_terms_map_.emplace( term_uuid, term_idx );
            return term_uuid;
        }

        [[nodiscard]] bool remove_energy_term( const uuid& term_id )
        {
            auto term_it = energy_terms_map_.find( term_id );
            if( term_it == energy_terms_map_.end() )
            {
                return false;
            }
            index_t idx = term_it->second;
            index_t last = energy_terms_.size() - 1;
            std::swap( energy_terms_[idx], energy_terms_[last] );
            energy_terms_map_[energy_terms_[idx]->id()] = idx;
            energy_terms_.pop_back();
            energy_terms_map_.erase( term_it );
            return true;
        }

        void clear()
        {
            energy_terms_.clear();
            energy_terms_map_.clear();
        }

        [[nodiscard]] index_t size() const
        {
            return energy_terms_.size();
        }

        [[nodiscard]] const EnergyTerm< ObjectType >& get(
            const uuid& term_id ) const
        {
            auto term_it = energy_terms_map_.find( term_id );
            OPENGEODE_EXCEPTION( term_it != energy_terms_map_.end(),
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    term_id.string() ) );
            return *energy_terms_[term_it->second];
        }

        [[nodiscard]] const std::vector<
            std::unique_ptr< EnergyTerm< ObjectType > > >&
            energy_terms() const
        {
            return energy_terms_;
        }

        [[nodiscard]] std::string string() const
        {
            auto message = absl::StrCat(
                "EnergyTermCollection: ", energy_terms_.size(), " terms:" );
            for( const auto& term : energy_terms_ )
            {
                absl::StrAppend( &message, "\n\t --> ", term->string() );
            }
            return message;
        }

    private:
        std::vector< std::unique_ptr< EnergyTerm< ObjectType > > >
            energy_terms_;
        absl::flat_hash_map< uuid, index_t > energy_terms_map_;
    };

} // namespace geode