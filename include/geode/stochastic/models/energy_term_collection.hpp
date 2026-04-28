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

        uuid add_energy_term(
            std::unique_ptr< EnergyTerm< ObjectType > >&& term )
        {
            auto term_idx = energy_terms_.size();

            const auto term_name = term->name();
            OPENGEODE_EXCEPTION( term_name.has_value(),
                absl::StrCat( "[EnergyTermCollection]- Energy Term name is not "
                              "defined." ) );
            const auto term_uuid = term->id();
            auto [it, inserted_uuid] =
                name_to_uuid_.emplace( term_name.value(), term_uuid );
            OPENGEODE_EXCEPTION( inserted_uuid,
                absl::StrCat( "[EnergyTermCollection]- Energy Term named ",
                    term_name.value(), " already exists." ) );

            auto [it2, inserted_index] =
                uuid_to_index_.emplace( term_uuid, term_idx );
            OPENGEODE_EXCEPTION( inserted_index,
                absl::StrCat( "[EnergyTermCollection]- Energy Term  ",
                    term_uuid.string(), " already exists." ) );

            energy_terms_.emplace_back( std::move( term ) );
            return term_uuid;
        }

        [[nodiscard]] index_t size() const
        {
            return energy_terms_.size();
        }

        [[nodiscard]] const EnergyTerm< ObjectType >& get(
            const uuid& term_id ) const
        {
            auto term_it = uuid_to_index_.find( term_id );
            OPENGEODE_EXCEPTION( term_it != uuid_to_index_.end(),
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    term_id.string() ) );
            return *energy_terms_[term_it->second];
        }

        [[nodiscard]] uuid get_term_uuid( std::string_view name ) const
        {
            auto uuid_it = name_to_uuid_.find( name );
            OPENGEODE_EXCEPTION( uuid_it != name_to_uuid_.end(),
                absl::StrCat(
                    "[EnergyTermCollection] Unknown energy term: ", name ) );
            return uuid_it->second;
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
        absl::flat_hash_map< std::string, uuid > name_to_uuid_;
        absl::flat_hash_map< uuid, index_t > uuid_to_index_;
        std::vector< std::unique_ptr< EnergyTerm< ObjectType > > >
            energy_terms_;
    };

} // namespace geode