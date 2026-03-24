#pragma once

#include <absl/container/flat_hash_map.h>
#include <geode/stochastic/sampling/mcmc/energy_terms/energy_term.hpp>

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
            std::shared_ptr< EnergyTerm< ObjectType > > term )
        {
            const uuid term_id = term->id();
            energy_terms_.emplace( term_id, term );
            for( const uuid& set_id : term->targeted_set_ids() )
            {
                set_to_terms_[set_id].push_back( term );
            }
            return term_id;
        }

        [[nodiscard]] bool remove_energy_term( const uuid& term_id )
        {
            auto term_it = energy_terms_.find( term_id );
            if( term_it == energy_terms_.end() )
            {
                return false;
            }

            auto term = term_it->second;

            for( const uuid& set_id : term->targeted_set_ids() )
            {
                auto vec_it = set_to_terms_.find( set_id );
                if( vec_it == set_to_terms_.end() )
                {
                    continue;
                }

                auto& vec = vec_it->second;
                vec.erase(
                    std::remove( vec.begin(), vec.end(), term ), vec.end() );

                if( vec.empty() )
                {
                    set_to_terms_.erase( vec_it );
                }
            }

            energy_terms_.erase( term_it );
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
            const uuid& term_id ) const
        {
            auto term_it = energy_terms_.find( term_id );
            OPENGEODE_EXCEPTION( term_it != energy_terms_.end(),
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    term_id.string() ) );
            return *term_it->second;
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

        [[nodiscard]] std::string string() const
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