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

add_geode_test(
    SOURCE "configuration/test-pairwise-interactions.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/direct/test-ball-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/direct/test-bounding-box-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/direct/test-double-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/direct/test-point-uniform-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/direct/test-segment-uniform-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/mcmc/models/components/test-intensity-term.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/mcmc/models/components/test-pairwise-term.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/mcmc/models/test-gibbs-energy.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/mcmc/proposal/test-proposal-kernel.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "sampling/mcmc/test-metropolis-hasting-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)
add_geode_test(
    SOURCE  "sampling/mcmc/test-mh-poisson-statistics-and-convergence.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
    UNSTABLE
)

add_geode_test(
    SOURCE  "sampling/mcmc/test-mh-strauss-statistics-and-convergence.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
    UNSTABLE
)


add_geode_test(
    SOURCE "sampling/test-random-engine.cpp"
    DEPENDENCIES
        OpenGeode::basic
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "test-hello-world.cpp"
    DEPENDENCIES
        OpenGeode::basic
        ${PROJECT_NAME}::stochastic
)



