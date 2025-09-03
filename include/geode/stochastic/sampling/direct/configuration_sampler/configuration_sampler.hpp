/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

namespace geode
{
    template < typename Object >
    class ConfigurationSampler
    {
    public:
        ConfigurationSampler( GroupId group_id ) : group_id_{ group_id } {}
        virtual ~ConfigurationSampler() = default;

        virtual Object sample( RandomEngine& engine ) const = 0;

        std::optional< ObjectId > sample_id(
            const Configuration< Object >& config, RandomEngine& engine ) const
        {
            const auto max_obj_id = config.nb_objects( group_id_ );
            if( max_obj_id == 0 )
            {
                return std::nullopt;
            }
            geode::UniformClosed< index_t > uniform_closed_index_t;
            uniform_closed_index_t.min_value = 0;
            uniform_closed_index_t.max_value = max_obj_id - 1;
            return { engine.sample_uniform( uniform_closed_index_t ),
                group_id_ };
        }

        virtual Object change(
            const Object& object, RandomEngine& engine ) const = 0;

        virtual double log_pdf( const MarkedObj& obj ) const = 0;

    private:
        GroupId group_id_;
    };
} // namespace geode