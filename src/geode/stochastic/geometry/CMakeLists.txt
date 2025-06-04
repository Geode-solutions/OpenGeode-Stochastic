#
# Copyright (c) 2019 - 2025 Geode-solutions. All rights reserved.
#

add_geode_library(
    NAME geometry
    FOLDER "geode/stochastic/geometry"
    SOURCES
        "common.cpp"
        "hello_world.cpp"
    PUBLIC_HEADERS
        "common.hpp"
        "hello_world.hpp"
    PRIVATE_DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry

)