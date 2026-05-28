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

#include <absl/container/btree_map.h>
#include <absl/strings/str_join.h>

#include <geode/basic/uuid.hpp>

#include <string>
#include <variant>
#include <vector>

#include <geode/stochastic/models/energy_term_collection.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_builder.hpp>
#include <geode/stochastic/models/energy_terms/energy_term_config.hpp>

#include <geode/stochastic/models/gibbs_energy.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>
namespace geode
{

    struct ModelConfig
    {
        std::vector< EnergyTermConfig > terms;
    };

    template < typename ObjectType >
    class Model
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( Model );

    public:
        Model() = delete;
        Model( EnergyTermCollection< ObjectType >&& energy_terms )
            : terms_collection_( std::move( energy_terms ) ),
              energy_{ terms_collection_ }
        {
        }

        [[nodiscard]] index_t nb_terms() const
        {
            return terms_collection_.size();
        }

        [[nodiscard]] const EnergyTermCollection< ObjectType >& terms() const
        {
            return terms_collection_;
        }

        [[nodiscard]] index_t term_index( const uuid& term_uuid ) const
        {
            return terms_collection_.get_term_index( term_uuid );
        }

        [[nodiscard]] const GibbsEnergy< ObjectType >& energy() const
        {
            return energy_;
        }

        [[nodiscard]] std::vector< double > compute_statistics(
            const ObjectSets< ObjectType >& state ) const
        {
            std::vector< double > statistic_values;
            statistic_values.reserve( terms_collection_.size() );
            for( const auto& term : terms_collection_.energy_terms() )
            {
                statistic_values.emplace_back( term->statistic( state ) );
            }
            return statistic_values;
        }

        [[nodiscard]] double compute_statistic(
            const ObjectSets< ObjectType >& state, const uuid& term_uuid ) const
        {
            const auto& term = terms_collection_.get( term_uuid );
            return term.statistic( state );
        }

        [[nodiscard]] std::string term_name( const uuid& term_uuid ) const
        {
            const auto& term = terms_collection_.get( term_uuid );
            return term.name().value_or( "unnamed" );
        }

        [[nodiscard]] std::vector< std::string > term_names() const
        {
            std::vector< std::string > names;
            names.reserve( terms_collection_.size() );

            for( const auto& term : terms_collection_.energy_terms() )
            {
                names.emplace_back( term->name().value_or( "unnamed" ) );
            }

            return names;
        }

        [[nodiscard]] std::string string() const
        {
            return terms_collection_.string();
        }

    private:
        EnergyTermCollection< ObjectType > terms_collection_;
        GibbsEnergy< ObjectType > energy_;
    };

    template < typename ObjectType >
    std::unique_ptr< Model< ObjectType > > build_model(
        const ModelConfig& config,
        const ObjectSets< ObjectType >& object_sets,
        const SpatialDomain< ObjectType::dim >& domain )
    {
        EnergyTermCollection< ObjectType > collection;

        for( const auto& term_cfg : config.terms )
        {
            collection.add_energy_term( build_energy_term< ObjectType >(
                term_cfg, object_sets, domain ) );
        }

        return std::make_unique< Model< ObjectType > >(
            std::move( collection ) );
    }
} // namespace geode