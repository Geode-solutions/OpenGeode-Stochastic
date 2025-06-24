#include <geode/stochastic/geometry/common.hpp>

#include <geode/stochastic/geometry/distributions.hpp>
#include <geode/stochastic/geometry/random_engine.hpp>

#include <geode/geometry/point.hpp>
#include <variant>

namespace geode
{
    using DoubleDistribution = std::variant< Uniform< double >, Gaussian >;
    using Domain = std::variant< UniformBox3D >;

    // Interface générique
    template < typename T >
    struct Sampler;

    template <>
    struct Sampler< double >
    {
        static double sample(
            RandomEngine& engine, const DoubleDistribution& dist )
        {
            return std::visit(
                [&engine]( auto&& d ) -> double {
                    using D = std::decay_t< decltype( d ) >;
                    if constexpr( std::is_same_v< D, Uniform< double > > )
                        return engine.sample_uniform< double >( d );
                    else if constexpr( std::is_same_v< D, Gaussian > )
                        return engine.sample_gaussian( d );
                    else
                        throw std::runtime_error(
                            "Unsupported distribution for double" );
                },
                dist );
        }
    };

    // Spécialisation pour `Point2D`
    template <>
    struct Sampler< geode::Point3D >
    {
        static geode::Point3D sample( RandomEngine& engine, const Domain& dist )
        {
            return std::visit(
                [&engine]( auto&& d ) -> geode::Point3D {
                    using D = std::decay_t< decltype( d ) >;
                    if constexpr( std::is_same_v< D, UniformBox3D > )
                        return engine.sample_box_uniform( d );
                    else
                        throw std::runtime_error(
                            "Unsupported distribution for Point2D" );
                },
                dist );
        }
    };
} // namespace geode