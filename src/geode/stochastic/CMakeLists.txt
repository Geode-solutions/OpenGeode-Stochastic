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
        "spatial/object_set.cpp"
        "spatial/object_sets.cpp"
        "spatial/object_neighborhood.cpp"
        "spatial/pairwise_interactions.cpp"
        "sampling/direct/ball_sampler.cpp"
        "sampling/direct/bounding_box_sampler.cpp"
        "sampling/direct/double_sampler.cpp"
        "sampling/direct/point_uniform_sampler.cpp"
        "sampling/direct/segment_uniform_sampler.cpp"
        "sampling/distributions.cpp"
        "sampling/random_engine.cpp"
        "sampling/mcmc/energy_terms/intensity_term.cpp"
        "common.cpp"
    PUBLIC_HEADERS
        "inference/abc_shadow.hpp"
        "models/spatial_model.hpp"
        "spatial/object_set.hpp"
        "spatial/object_sets.hpp"
        "spatial/object_neighborhood.hpp"
        "spatial/object_helpers.hpp"
        "spatial/pairwise_interactions.hpp"
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
        "sampling/mcmc/helpers/fracture_simulation_runner.hpp"
        "sampling/mcmc/helpers/simulation_printer.hpp"
        "sampling/mcmc/helpers/simulation_monitor.hpp"
        "sampling/mcmc/energy_terms/density_term.hpp"
        "sampling/mcmc/energy_terms/energy_term.hpp"
        "sampling/mcmc/energy_terms/intensity_term.hpp"
        "sampling/mcmc/energy_terms/pairwise_term.hpp"
        "sampling/mcmc/energy_terms/single_object_term.hpp"
        "sampling/mcmc/energy_terms/energy_term_collection.hpp"
        "sampling/mcmc/energy_terms/gibbs_energy.hpp"
        "sampling/mcmc/proposal/classical_proposals.hpp"
        "sampling/mcmc/proposal/moves.hpp"
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