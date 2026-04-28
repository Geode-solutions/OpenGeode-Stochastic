#pragma once
#include <geode/basic/common.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/stochastic/common.hpp>

namespace geode
{

    class StatMonitor
    {
    public:
        void add_realization(
            const absl::flat_hash_map< uuid, double >& values )
        {
            ++count_;
            for( const auto& [term_uuid, value] : values )
            {
                auto& mean = means_[term_uuid];
                auto& m2 = m2_[term_uuid]; // somme des carrés

                double delta = value - mean;
                mean += delta / count_;
                double delta2 = value - mean;
                m2 += delta * delta2;
            }
        }

        const index_t statiscal_count() const
        {
            return count_;
        }

        double mean( const uuid& id ) const
        {
            return means_.at( id );
        }

        double variance( const uuid& id ) const
        {
            if( count_ < 2 )
                return 0.0;
            return m2_.at( id ) / ( count_ - 1 );
        }

    private:
        absl::flat_hash_map< uuid, double > means_;
        absl::flat_hash_map< uuid, double > m2_;
        index_t count_{ 0 };
    };
} // namespace geode