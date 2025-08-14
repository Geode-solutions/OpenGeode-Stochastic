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

#include <geode/stochastic/common.hpp>

#include <geode/stochastic/configuration/configuration.hpp>

namespace
{
    struct NegLogParam
    {
        explicit NegLogParam( double param )
        {
            OPENGEODE_EXCEPTION( param >= 0.,
                "[Gibbs energy term] - The model parameter "
                "cannot be smaller than 0." );
            if( param >= geode::GLOBAL_EPSILON )
            {
                value = -std::log( param );
            }
        }

        double scale( double multiplier ) const
        {
            if( value )
            {
                return value.value() * multiplier;
            }
            return ( multiplier > 0 )
                       ? std::numeric_limits< double >::infinity()
                       : 0.0;
        }

        double param() const
        {
            if( value )
            {
                return std::exp( -value.value() );
            }
            return 0.;
        }
        std::optional< double > value;
    };
} // namespace

namespace geode
{
    template < typename Geometry >
    class EnergyTerm
    {
    public:
        explicit EnergyTerm( double parameter )
            : neg_log_parameter_( parameter )
        {
        }

        virtual ~EnergyTerm() = default;

        double parameter() const
        {
            return neg_log_parameter_.param();
        }

        virtual double log_total(
            const Configuration< Geometry >& state ) const = 0;

        virtual double log_delta_add( const Configuration< Geometry >& state,
            const MarkedObject< Geometry >& sample ) const = 0;

        virtual double log_delta_remove( const Configuration< Geometry >& state,
            index_t sample_id ) const = 0;

        virtual double log_delta_change( const Configuration< Geometry >& state,
            index_t old_sample_id,
            const MarkedObject< Geometry >& new_sample ) const = 0;

        virtual double statistic(
            const Configuration< Geometry >& state ) const = 0;

    protected:
        NegLogParam neg_log_parameter_;
    };
} // namespace geode