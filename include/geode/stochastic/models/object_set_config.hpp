
class ObjectSetModelConfig
{
public:
    virtual ~ObjectSetModelConfig() = default;
    virtual void add_energy_terms( const uuid& set_id,
        EnergyTermCollection< OwnerSegment2D >& collection,
        const SpatialDomain< 2 >& domain ) const = 0;
};

class PoissonSetModelConfig : public ObjectSetModelConfig
{
public:
    double intensity{ 1.0 };
    virtual void add_energy_terms( const uuid& set_id,
        EnergyTermCollection< OwnerSegment2D >& collection,
        const SpatialDomain< 2 >& domain ) const
    {
        collection.add_energy_term(
            std::make_unique< geode::DensityTerm< geode::Point2D > >(
                absl::StrCat( energy_desc.name, "_density" ),
                energy_desc.density, std::vector< geode::uuid >{ set_id },
                this->domain_ ) );
    }
};

class ObjectSetDynamicsConfig
{
public:
    double birth_weight{ 1.0 };
    double death_weight{ 1.0 };
    double change_weight{ 1.0 };
};

class ObjectSetConfig
{
public:
    std::string name;

    std::unique_ptr< ObjectSetModelConfig > set_model_;
    SetDynamicsConfig dynamics;
    std::unique_ptr< ObjectSetSampler< OwnerSegment2D > > build_sampler(
        const SpatialDomain< 2 >& domain ) const
    {
        return std::make_unique< UniformSegmentSetSampler >( domain );
    }
};

class ObjectSetsConfig
{
public:
    std::string name;
    // Nom du type de modèle
    virtual std::string model_name() const = 0;
    std::vector< ObjectSetConfig > objectset_config_;

    EnergyTermCollection< OwnerSegment2D > build_energy_terms(
        const SpatialDomain< 2 >& domain ) const
    {
        EnergyTermCollection< OwnerSegment2D > collection;
        for( const auto& set : sets )
        {
            const uuid set_id = uuid::create();
            set.model->add_energy_terms( set_id, collection, domain );
        }
        return collection;
    }

    std::vector< std::unique_ptr< ObjectSetSampler< OwnerSegment2D > > >
        build_samplers( const SpatialDomain< 2 >& domain ) const
    {
        std::vector< std::unique_ptr< ObjectSetSampler< OwnerSegment2D > > >
            samplers;
        for( const auto& set : sets )
        {
            samplers.push_back( set.dynamics.build_sampler( domain ) );
        }
        return samplers;
    }

    std::unique_ptr< ProposalKernel< OwnerSegment2D > > build_kernel() const
    {
        auto kernel = std::make_unique< ProposalKernel< OwnerSegment2D > >();
        for( const auto& set : sets )
        {
            const uuid set_id = uuid::create(); // doit correspondre au set réel
            kernel->add_move(
                set_id, std::make_unique< BirthMove< OwnerSegment2D > >(
                            set.dynamics.birth_weight ) );
            kernel->add_move(
                set_id, std::make_unique< DeathMove< OwnerSegment2D > >(
                            set.dynamics.death_weight ) );
            kernel->add_move(
                set_id, std::make_unique< ChangeMove< OwnerSegment2D > >(
                            set.dynamics.change_weight ) );
        }
        return kernel;
    }
};
class FractureModelConfig : public ObjectSetsConfig
{
public:
    double beta_x_node_{ 1. };
};