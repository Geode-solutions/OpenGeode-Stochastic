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

add_subdirectory(stochastic)


# The stubs of the binding module and of the pure-Python module must pass the
# quality gates (named parameters, no unresolved or raw C++ types, stubtest).
# The stubs of the OpenGeode modules whose types appear in the signatures are
# generated with it, for mypy to resolve them
add_geode_python_stub_test(
    MODULES
        opengeode_py_basic
        opengeode_py_geometry
        opengeode_stochastic_py_stochastic
    DEPENDENCIES
        ${PROJECT_NAME}::py_stochastic
)

add_geode_python_stub_test(
    NAME test-py-wrapper-stubs
    NO_STUBTEST
    MODULES
        stochastic
    DEPENDENCIES
        ${PROJECT_NAME}::py_stochastic
)
