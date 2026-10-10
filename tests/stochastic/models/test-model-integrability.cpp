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

namespace
{
    geode::PairwiseTermConfig pairwise_term( std::string_view name,
        std::vector< std::pair< std::string, std::string > > set_names,
        double gamma,
        double distance )
    {
        geode::PairwiseTermConfig term;
        term.term_name = name;
        term.object_set_names_interactions = std::move( set_names );
        term.gamma = gamma;
        term.interaction_config =
            geode::MinimalDistanceCutoffConfig{ distance };
        return term;
    }

    bool is_integrable( const geode::ModelConfig& config )
    {
        try
        {
            geode::check_model_integrability( config );
            return true;
        }
        catch( const geode::OpenGeodeException& )
        {
            return false;
        }
    }

    void check( const geode::ModelConfig& config,
        bool expected_integrable,
        std::string_view message )
    {
        geode::OpenGeodeStochasticStochasticException::test(
            is_integrable( config ) == expected_integrable,
            "[ModelIntegrability] ", message );
    }

    void test_repulsive_model()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::ModelConfig config;
        config.terms.emplace_back(
            pairwise_term( "repulsion", { { "A", "A" } }, 0.5, 1. ) );
        config.terms.emplace_back(
            pairwise_term( "x_node", { { "A", "B" } }, 0., 0. ) );
        // NOLINTEND(*-magic-numbers)
        check( config, true, "repulsive model should be integrable" );
    }

    void test_attractive_model_without_hard_core()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::ModelConfig config;
        config.terms.emplace_back(
            pairwise_term( "attraction", { { "A", "A" } }, 2., 1. ) );
        // NOLINTEND(*-magic-numbers)
        check( config, false,
            "attractive model without hard-core should not be integrable" );
    }

    void test_attractive_model_with_zero_hard_core()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::ModelConfig config;
        config.terms.emplace_back(
            pairwise_term( "attraction", { { "A", "A" } }, 2., 1. ) );
        // non-intersection does not bound the number of objects
        config.terms.emplace_back(
            pairwise_term( "spacing", { { "A", "A" } }, 0., 0. ) );
        // NOLINTEND(*-magic-numbers)
        check( config, false,
            "attractive model with a zero distance hard-core should not be "
            "integrable" );
    }

    void test_attractive_model_with_hard_core()
    {
        // NOLINTBEGIN(*-magic-numbers)
        geode::ModelConfig config;
        config.terms.emplace_back(
            pairwise_term( "attraction", { { "A", "A" } }, 2., 1. ) );
        config.terms.emplace_back(
            pairwise_term( "spacing", { { "A", "A" } }, 0., 0.1 ) );
        // NOLINTEND(*-magic-numbers)
        check( config, true,
            "attractive model with hard-core should be integrable" );
    }

    void test_attractive_inter_set_model()
    {
        // NOLINTBEGIN(*-magic-numbers)
        // e.g. favored Y-nodes between two fracture sets
        geode::ModelConfig config;
        config.terms.emplace_back(
            pairwise_term( "y_node", { { "A", "B" } }, 3., 0. ) );
        config.terms.emplace_back(
            pairwise_term( "spacing_A", { { "A", "A" } }, 0., 0.5 ) );
        check( config, false,
            "inter-set attraction with a single hard-core set should not be "
            "integrable" );

        config.terms.emplace_back(
            pairwise_term( "spacing_B", { { "B", "B" } }, 0., 0.5 ) );
        check( config, true,
            "inter-set attraction with hard-core on both sets should be "
            "integrable" );

        // an inter-set hard-core does not bound the number of objects per set
        geode::ModelConfig inter_hard_core;
        inter_hard_core.terms.emplace_back(
            pairwise_term( "y_node", { { "A", "B" } }, 3., 0. ) );
        inter_hard_core.terms.emplace_back(
            pairwise_term( "spacing_AB", { { "A", "B" } }, 0., 0.5 ) );
        check( inter_hard_core, false,
            "inter-set hard-core should not make the model integrable" );
        // NOLINTEND(*-magic-numbers)
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeStochasticStochasticLibrary::initialize();
        geode::Logger::info( "TEST - MODEL INTEGRABILITY" );

        test_repulsive_model();
        test_attractive_model_without_hard_core();
        test_attractive_model_with_zero_hard_core();
        test_attractive_model_with_hard_core();
        test_attractive_inter_set_model();

        geode::Logger::info( "--> SUCCESS!" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
