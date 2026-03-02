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

namespace geode
{
    template < typename Type >
    class ObjectSetSampler
    {
    public:
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
        virtual bool is_valid_object( const Type& obj ) const = 0;

    protected:
        double log_pdf_{ LOG_PROB_INVALID };
    };
} // namespace geode