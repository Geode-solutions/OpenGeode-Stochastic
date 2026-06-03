#pragma once

#include <geode/stochastic/inference/target_statistics.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/segment_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_context.hpp>
#include <geode/stochastic/sampling/mcmc/simulation_runner.hpp>

namespace geode
{
    using Fracture = OwnerSegment2D;
    using FractureSamplerConfig = ObjectSamplerConfig< Fracture >;
    using FractureSimulationContext = SimulationContext< Fracture >;
    using FractureSimulationRunner = SimulationRunner< Fracture >;

    struct FractureSetDescription
    {
        std::string fset_name;

        FractureSamplerConfig sampler;
        double birth_ratio{ 1.0 };
        double death_ratio{ 1.0 };
        double change_ratio{ 1.0 };

        std::string density_name() const
        {
            return absl::StrCat( fset_name, "_p20" );
        }
        double p20{ 0. };
        std::optional< double > expected_number;

        std::string intensity_name() const
        {
            return absl::StrCat( fset_name, "_p21" );
        }
        double p21{ 0. };
        std::optional< double > expected_total_length;

        std::vector< std::array< geode::Point2D, 2 > > observed_fractures;

        std::string spacing_name() const
        {
            return absl::StrCat( fset_name, "_spacing" );
        }
        double minimal_spacing{ 0. };

        std::string string() const
        {
            auto message =
                absl::StrCat( "FractureSetDescription: ", fset_name );
            for( const auto& fixed_object : observed_fractures )
            {
                absl::StrAppend( &message,
                    "\n\t --> observation (x,y,z)start: ",
                    fixed_object[0].string(),
                    " (x,y,z)end: ", fixed_object[1].string() );
            }
            absl::StrAppend( &message,
                "\n\t --> length distribution: ", sampler.length.string() );
            absl::StrAppend( &message,
                "\n\t --> azimuth distribution: ", sampler.azimuth.string() );
            absl::StrAppend( &message, "\n\t --> ", density_name(), ": ", p20 );
            absl::StrAppend(
                &message, "\n\t --> ", intensity_name(), ": ", p21 );
            absl::StrAppend(
                &message, "\n\t --> ", spacing_name(), ": ", minimal_spacing );

            absl::StrAppend( &message,
                "\n\t --> dynamic move ratio - birth/death/change (",
                birth_ratio, " / ", death_ratio, " / ", change_ratio, ")" );
            return message;
        }
    };

    //    struct FractureInterSetDescription
    //    {
    //        std::string interaction_name;
    //
    //        std::vector< std::string > set_names;
    //
    //        double gamma{ 1. };
    //        double distance{ 0. };
    //
    //        bool include_intra_set{ true };
    //        bool include_inter_set{ false };
    //
    //        std::optional< double > expected_nb_interactions;
    //    };

    struct FractureNetworkDescription
    {
        std::string fnet_name;

        SpatialDomainConfig< 2 > domain;

        std::vector< FractureSetDescription > fracture_sets;

        FractureSetDescription& add_fracture_set( absl::string_view fset_name )
        {
            auto& fracture_set = fracture_sets.emplace_back();
            fracture_set.fset_name = fset_name;
            return fracture_set;
        }

        std::string string() const
        {
            auto message =
                absl::StrCat( "FractureNetworkDescription: ", fnet_name );
            absl::StrAppend( &message, "\n\t --> ", domain.string() );
            for( const auto& fset_desc : fracture_sets )
            {
                absl::StrAppend( &message, "\n\t --> ", fset_desc.string() );
            }
            return message;
        }
        //        std::vector< StraussInteractionDescription< ObjectType > >
        //            inter_set_interactions;
        //
        //        void add_x_node_monitoring( double beta_x_node )
        //        {
        //            OpenGeodeStochasticStochasticException::check_exception(
        //                beta_x_node <= 1.0 && beta_x_node >= 0., nullptr,
        //                OpenGeodeException::TYPE::data,
        //                "[FractureSimulationRunner] x node should be
        //                inhibitated, " "please provise a value in [0., 1.]."
        //                );
        //            beta_x_node_ = beta_x_node;
        //        }
    };

    FractureSimulationContext build_fractures_simulation_context(
        const FractureNetworkDescription& description );

    std::vector< geode::TargetStatisticConfig > build_fractures_targeted_stat(
        const FractureNetworkDescription& description );

} // namespace geode