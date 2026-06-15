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

#include <geode/stochastic/sampling/random_engine.hpp>

#include <geode/stochastic/spatial/object_sets.hpp>
#include <geode/stochastic/spatial/spatial_domain.hpp>

namespace geode
{

    template < typename Type >
    class ObjectSetSampler
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( ObjectSetSampler );

    public:
        ObjectSetSampler() = default;
        virtual ~ObjectSetSampler() = default;

        [[nodiscard]] virtual Type sample( RandomEngine& engine ) const = 0;

        [[nodiscard]] virtual Type change(
            const Type& object, RandomEngine& engine ) const = 0;

        [[nodiscard]] double log_pdf( const Type& obj ) const
        {
            if( !is_valid_object( obj ) )
            {
                Logger::warn( "[ObjectSetSampler] - invalid object proposed." );
                return LOG_PROB_INVALID;
            }
            return log_pdf_;
        }

    protected:
        [[nodiscard]] virtual bool is_valid_object( const Type& obj ) const = 0;
        void set_log_pdf( double value )
        {
            log_pdf_ = value;
        }

    private:
        double log_pdf_{ LOG_PROB_INVALID };
    };

    template < typename ObjectType >
    struct ObjectSamplerConfig;

    template < typename ObjectType >
    std::unique_ptr< ObjectSetSampler< ObjectType > > build_objectset_sampler(
        const SpatialDomain< ObjectType::dim >& domain,
        const ObjectSamplerConfig< ObjectType >& config );
} // namespace geode