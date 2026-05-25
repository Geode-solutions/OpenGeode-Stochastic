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
            OpenGeodeStochasticStochasticException::check_exception(
                term_name.has_value(), nullptr, OpenGeodeException::TYPE::data,
                absl::StrCat( "[EnergyTermCollection]- Energy Term name is not "
                              "defined." ) );
            const auto term_uuid = term->id();
            auto [it, inserted_uuid] =
                name_to_uuid_.emplace( term_name.value(), term_uuid );
            OpenGeodeStochasticStochasticException::check_exception(
                inserted_uuid, nullptr, OpenGeodeException::TYPE::data,
                absl::StrCat( "[EnergyTermCollection]- Energy Term named ",
                    term_name.value(), " already exists." ) );

            auto [it2, inserted_index] =
                uuid_to_index_.emplace( term_uuid, term_idx );
            OpenGeodeStochasticStochasticException::check_exception(
                inserted_index, nullptr, OpenGeodeException::TYPE::data,
                absl::StrCat( "[EnergyTermCollection]- Energy Term  ",
                    term_uuid.string(), " already exists." ) );

            energy_terms_.emplace_back( std::move( term ) );
            return term_uuid;
        }

        [[nodiscard]] index_t size() const
        {
            return energy_terms_.size();
        }

        [[nodiscard]] index_t get_term_index( const uuid& term_uuid ) const
        {
            auto term_it = uuid_to_index_.find( term_uuid );
            OpenGeodeStochasticStochasticException::check_exception(
                term_it != uuid_to_index_.end(), nullptr,
                OpenGeodeException::TYPE::data,
                absl::StrCat( "[EnergyTermCollection] Unknown energy term: ",
                    term_uuid.string() ) );
            return term_it->second;
        }

        [[nodiscard]] const EnergyTerm< ObjectType >& get(
            const uuid& term_uuid ) const
        {
            return *energy_terms_[get_term_index( term_uuid )];
        }

        [[nodiscard]] uuid get_term_uuid( std::string_view name ) const
        {
            auto uuid_it = name_to_uuid_.find( name );
            OpenGeodeStochasticStochasticException::check_exception(
                uuid_it != name_to_uuid_.end(), nullptr,
                OpenGeodeException::TYPE::data,
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