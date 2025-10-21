# Copyright (c) 2019 - 2025 Geode-solutions
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
        "sampling/direct/ball_sampler.cpp"
        "sampling/direct/bounding_box_sampler.cpp"
        "sampling/direct/double_sampler.cpp"
        "sampling/direct/point_uniform_sampler.cpp"
        "sampling/direct/segment_uniform_sampler.cpp"
        "sampling/distributions.cpp"
        "sampling/random_engine.cpp"
        "common.cpp"
        "hello_world.cpp"
    PUBLIC_HEADERS
        "spatial/object_set.hpp"
        "spatial/object_sets.hpp"
        "spatial/object_neighborhood.hpp"
        "spatial/object_helpers.hpp"
        "spatial/pairwise_interactions.hpp"
        "spatial/details/RTree.hpp"
        "sampling/direct/object_set_sampler/object_set_sampler.hpp"
        "sampling/direct/object_set_sampler/point_set_sampler.hpp"
        "sampling/direct/ball_sampler.hpp"
        "sampling/direct/bounding_box_sampler.hpp"
        "sampling/direct/double_sampler.hpp"
        "sampling/direct/point_uniform_sampler.hpp"
        "sampling/direct/segment_uniform_sampler.hpp"
        "sampling/mcmc/helpers/simulation_runner.hpp"
        "sampling/mcmc/models/components/pairwise_term.hpp"
        "sampling/mcmc/models/components/density_term.hpp"
        "sampling/mcmc/models/components/energy_term.hpp"
        "sampling/mcmc/models/energy_term_collection.hpp"
        "sampling/mcmc/models/gibbs_energy.hpp"
        "sampling/mcmc/proposal/classical_proposals.hpp"
        "sampling/mcmc/proposal/moves.hpp"
        "sampling/mcmc/proposal/proposal_kernel.hpp"
        "sampling/mcmc/metropolis_hasting_sampler.hpp"
        "sampling/distributions.hpp"
        "sampling/random_engine.hpp"
        "common.hpp"
        "hello_world.hpp"
    PRIVATE_DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        absl::random_random
        absl::random_distributions

)