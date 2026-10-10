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

#include <geode/stochastic/models/model.hpp>

#include <absl/container/flat_hash_set.h>

namespace
{
    std::optional< double > interaction_distance(
        const geode::PairwiseInteractionConfig& interaction_config )
    {
        return std::visit(
            []( const auto& config ) -> std::optional< double > {
                using Config = std::decay_t< decltype( config ) >;
                if constexpr( std::is_same_v< Config, std::monostate > )
                {
                    return std::nullopt;
                }
                else
                {
                    return config.threshold;
                }
            },
            interaction_config );
    }

    bool is_hard_core( const geode::PairwiseTermConfig& term )
    {
        // gamma below epsilon is handled as 0 by the energy scale
        if( term.gamma >= geode::GLOBAL_EPSILON )
        {
            return false;
        }
        const auto distance = interaction_distance( term.interaction_config );
        return distance && *distance > 0.;
    }

    // Sets with an intra-set hard-core at a strictly positive distance: their
    // number of objects is bounded in the (bounded) simulation domain.
    absl::flat_hash_set< std::string > hard_core_sets(
        const geode::ModelConfig& config )
    {
        absl::flat_hash_set< std::string > sets;
        for( const auto& term_cfg : config.terms )
        {
            const auto* term =
                std::get_if< geode::PairwiseTermConfig >( &term_cfg );
            if( term == nullptr || !is_hard_core( *term ) )
            {
                continue;
            }
            for( const auto& [set1, set2] :
                term->object_set_names_interactions )
            {
                if( set1 == set2 )
                {
                    sets.insert( set1 );
                }
            }
        }
        return sets;
    }
} // namespace

namespace geode
{
    void check_model_integrability( const ModelConfig& config )
    {
        const auto hard_cores = hard_core_sets( config );
        for( const auto& term_cfg : config.terms )
        {
            const auto* term = std::get_if< PairwiseTermConfig >( &term_cfg );
            if( term == nullptr || term->gamma <= 1. )
            {
                continue;
            }
            for( const auto& set_names : term->object_set_names_interactions )
            {
                for( const auto& set_name :
                    { set_names.first, set_names.second } )
                {
                    OpenGeodeStochasticStochasticException::check_exception(
                        hard_cores.contains( set_name ), nullptr,
                        OpenGeodeException::TYPE::data,
                        "[Model] Attractive pairwise term '", term->term_name,
                        "' (gamma = ", term->gamma,
                        " > 1) is not integrable: object set '", set_name,
                        "' requires a hard-core interaction (gamma = 0) at a "
                        "strictly positive distance." );
                }
            }
        }
    }
} // namespace geode
