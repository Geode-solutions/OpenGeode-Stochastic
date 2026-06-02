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

add_geode_test(
    SOURCE "direct/test-ball-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "direct/test-bounding-box-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "direct/test-double-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "direct/test-point-uniform-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

add_geode_test(
    SOURCE "direct/test-segment-uniform-sampler.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

#add_geode_test(
#    SOURCE "mcmc/helpers/test-simulation_printer.cpp"
#    DEPENDENCIES
#        OpenGeode::basic
#        OpenGeode::geometry
#        ${PROJECT_NAME}::stochastic
#)

#add_geode_test(
#    SOURCE "mcmc/helpers/test-simulation_monitor.cpp"
#    DEPENDENCIES
#        OpenGeode::basic
#        OpenGeode::geometry
#        ${PROJECT_NAME}::stochastic
#)

add_geode_test(
    SOURCE "mcmc/proposal/test-proposal-kernel.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        ${PROJECT_NAME}::stochastic
)

#add_geode_test(
#    SOURCE "mcmc/test-metropolis-hasting-sampler.cpp"
#    DEPENDENCIES
#        OpenGeode::basic
#        OpenGeode::geometry
#        ${PROJECT_NAME}::stochastic
#)

add_geode_test(
    SOURCE "test-random-engine.cpp"
    DEPENDENCIES
        OpenGeode::basic
        ${PROJECT_NAME}::stochastic
)
