#pragma once

// #include <absl/container/flat_hash_map.h>
#include <geode/basic/common.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/stochastic/common.hpp>
#include <geode/stochastic/models/model.hpp>

namespace geode
{
    template < typename ObjectType >
    class StatisticsTracker
    {
    public:
        StatisticsTracker( const Model< ObjectType >& model ) : model_{ model }
        {
            means_.resize( model.nb_terms(), 0.0 );
            m2_.resize( model.nb_terms(), 0.0 );
        }

        [[nodiscard]] index_t statiscal_count() const
        {
            return count_;
        }

        void add_realization( const std::vector< double >& values )
        {
            ++count_;
            for( const auto value_id : geode::Range{ values.size() } )
            {
                auto& value = values[value_id];
                auto& mean = means_[value_id];
                auto& m2 = m2_[value_id]; // somme des carrés

                double delta = value - mean;
                mean += delta / count_;
                double delta2 = value - mean;
                m2 += delta * delta2;
            }
        }

        [[nodiscard]] double mean( const uuid& term_uuid ) const
        {
            return means_[model_.term_index( term_uuid )];
        }

        [[nodiscard]] const std::vector< double >& means() const
        {
            return means_;
        }

        [[nodiscard]] double variance( const uuid& term_uuid ) const
        {
            return variance( model_.term_index( term_uuid ) );
        }

        [[nodiscard]] std::vector< double > variances() const
        {
            std::vector< double > variances;
            variances.reserve( model_.nb_terms() );
            for( const auto variance_id : geode::Range{ model_.nb_terms() } )
            {
                variances.emplace_back( this->variance( variance_id ) );
            }
            return variances;
        }

    private:
        [[nodiscard]] double variance( index_t term_index ) const
        {
            if( count_ < 2 )
            {
                return 0.0;
            }
            return m2_[term_index] / ( count_ - 1 );
        }

    private:
        const Model< ObjectType >& model_;

        std::vector< double > means_{};
        std::vector< double > m2_{};
        index_t count_{ 0 };
    };
} // namespace geode