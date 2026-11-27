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

#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>

namespace geode
{
    /// Meaning of the density given to an object set with fixed (observed)
    /// objects. For fractures, see FractureProcessBuilder::add_fracture_set
    /// for the choice according to what was observed.
    enum class IntensityScope
    {
        /// The density is the one of the simulated objects only: all the fixed
        /// objects are added to them.
        free_objects,
        /// The density is the one of all the objects of the length law of the
        /// set, fixed ones included. The fixed objects within the length law
        /// that touch the domain are taken as a random sample of them,
        /// independent of size and position: the simulated objects touching
        /// the domain get the remaining density lambda * (1 - p), with
        /// p = F / (lambda * A) (F such fixed objects, A from
        /// touching_domain_area), the others keep lambda. The other fixed
        /// objects (outside the length law, e.g. all the objects longer than
        /// an observation cutoff, or not touching the domain) are added as
        /// with free_objects. Exact for a Poisson process, an approximation
        /// with interaction terms on the set.
        all_objects
    };

    /// Area of the anchor positions for which an object drawn by the sampler
    /// touches the domain (SpatialDomainChecker::touches_domain).
    /// - Points: the area of the domain.
    /// - Segments: the mean area of the Minkowski sum of the domain and of the
    /// segments (Chiu et al., 2013), for a w x h domain
    /// w * h + E[L * (w |cos theta| + h |sin theta|)]. The expectation is
    /// estimated by Monte Carlo with a fixed seed on the length and azimuth
    /// laws only (any law): w * h is exact, and the relative error on the
    /// boundary term, about 0.2 %, does not depend on the domain size.
    template < typename ObjectType >
    [[nodiscard]] double touching_domain_area(
        const SpatialDomain< ObjectType::dim >& domain,
        const ObjectSamplerConfig< ObjectType >& sampler_config );

    /// Adds to the model the observation term of the set at set_index in
    /// config.sets, whose density lambda
    /// counts all the objects (IntensityScope::all_objects): a single object
    /// term with parameter 1 - p on the objects touching the domain
    /// (ObjectTouchingDomainFeature). Throws if lambda * A is lower than the
    /// number of fixed objects within the length law touching the domain. If
    /// they are equal, the free objects cannot touch the domain.
    template < typename ObjectType >
    void add_observation_term( SimulationContextConfig< ObjectType >& config,
        index_t set_index,
        double lambda );
} // namespace geode
