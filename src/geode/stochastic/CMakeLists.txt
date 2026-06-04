# Copyright (c) 2019 - 2026 Geode-solutions
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

add_geode_library(
    NAME stochastic
    FOLDER "geode/stochastic"
    SOURCES
        "applications/fractures.cpp"
        "applications/poisson_process.cpp"
        "applications/strauss_process.cpp"
        "spatial/object_set.cpp"
        "spatial/object_sets.cpp"
        "spatial/object_neighborhood.cpp"
        "spatial/single_object_features/segment_length_feature.cpp"
        "spatial/pairwise_interactions/distance_cutoff.cpp"
        "sampling/direct/object_set_sampler/segment_set_sampler.cpp"
        "sampling/direct/object_set_sampler/point_set_sampler.cpp"
        "sampling/direct/ball_sampler.cpp"
        "sampling/direct/bounding_box_sampler.cpp"
        "sampling/direct/double_sampler.cpp"
        "sampling/direct/point_uniform_sampler.cpp"
        "sampling/direct/segment_uniform_sampler.cpp"
        "models/energy_terms/energy_term_builder.cpp"
        "sampling/distributions.cpp"
        "sampling/random_engine.cpp"
        "common.cpp"
    PUBLIC_HEADERS
        "applications/fractures.hpp"
        "applications/poisson_process.hpp"
        "applications/strauss_process.hpp"
        #"inference/abc_shadow.hpp"
        "inference/statistics_tools.hpp"
        "inference/statistics_tracker.hpp"
        "inference/target_statistics.hpp"
        "spatial/object_set.hpp"
        "spatial/object_sets.hpp"
        "spatial/object_neighborhood.hpp"
        "spatial/object_helpers.hpp"
        "spatial/single_object_features/single_object_feature_builder.hpp"
        "spatial/single_object_features/single_object_feature_config.hpp"
        "spatial/single_object_features/segment_length_feature.hpp"
        "spatial/single_object_features/single_object_feature.hpp"
        "spatial/pairwise_interactions/pairwise_interactions_builder.hpp"
        "spatial/pairwise_interactions/pairwise_interactions_config.hpp"
        "spatial/pairwise_interactions/pairwise_interactions.hpp"
        "spatial/pairwise_interactions/distance_cutoff.hpp"
        "spatial/spatial_domain.hpp"
        "spatial/details/RTree.hpp"
        "sampling/direct/object_set_sampler/object_set_sampler.hpp"
        "sampling/direct/object_set_sampler/point_set_sampler.hpp"
        "sampling/direct/object_set_sampler/segment_set_sampler.hpp"
        "sampling/direct/ball_sampler.hpp"
        "sampling/direct/bounding_box_sampler.hpp"
        "sampling/direct/double_sampler.hpp"
        "sampling/direct/point_uniform_sampler.hpp"
        "sampling/direct/segment_uniform_sampler.hpp"
        #"sampling/mcmc/helpers/fracture_simulation_runner.hpp"
        "sampling/mcmc/helpers/simulation_context.hpp"
        "sampling/mcmc/helpers/simulation_printer.hpp"
        "models/energy_terms/energy_term.hpp"
        "models/energy_terms/pairwise_term.hpp"
        "models/energy_terms/single_object_term.hpp"
        "models/energy_terms/energy_term_config.hpp"
        "models/energy_terms/energy_term_builder.hpp"
        "models/energy_term_collection.hpp"
        "models/gibbs_energy.hpp"
        "models/model.hpp"
        "sampling/mcmc/proposal/classical_proposals.hpp"
        "sampling/mcmc/proposal/moves.hpp"
        "sampling/mcmc/proposal/object_set_dynamic_config.hpp"
        "sampling/mcmc/proposal/proposal_kernel.hpp"
        "sampling/mcmc/metropolis_hasting_sampler.hpp"
        "sampling/mcmc/simulation_runner.hpp"
        "sampling/distributions.hpp"
        "sampling/random_engine.hpp"
        "common.hpp"
    PRIVATE_DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        absl::random_random
        absl::random_distributions

)