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
    /// objects. Three observation cases, with lambda the density of the
    /// objects of the length law of the set:
    /// 1. a sample of the objects of all sizes: the length law of all the
    /// objects, all_objects;
    /// 2. all the objects longer than a cutoff L_obs, none shorter: a length
    /// law cut at L_obs, all_objects or free_objects (identical here);
    /// 3. all the objects longer than L_obs and a sample of the shorter ones:
    /// a length law cut at L_obs, all_objects.
    /// In cases 2 and 3, the fixed objects longer than L_obs act on the
    /// simulated ones through the interaction terms of the set.
    enum class IntensityScope
    {
        /// Density of the free objects only: the fixed objects are added to
        /// them, without changing their law.
        free_objects,
        /// Density of all the objects, fixed ones included. The fixed objects
        /// within the length law are a part of the objects touching the
        /// domain, observed independently of their size and position: an
        /// independent thinning of the Poisson process (Chiu et al., 2013).
        /// The objects which were not observed are then a Poisson process of
        /// density lambda * (1 - p) on the objects touching the domain, and
        /// lambda on the others, with p = F / (lambda * A): F fixed objects
        /// within the length law touch the domain, and lambda * A objects
        /// are expected to touch it (see touching_domain_area). Fixed objects
        /// outside the length law (longer than an observation cutoff, all
        /// observed) or not touching the domain are added to the free ones,
        /// as with free_objects. This is exact for a Poisson process; with
        /// interaction terms on the set, it is only an approximation, whose
        /// parameters should rather be inferred.
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
