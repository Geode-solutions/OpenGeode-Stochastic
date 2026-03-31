template < typename ObjectType >
std::unique_ptr< PairwiseInteraction< ObjectType > > build_interaction(
    const PairwiseInteractionConfig& cfg )
{
    return std::visit(
        [&](
            auto&& c ) -> std::unique_ptr< PairwiseInteraction< ObjectType > > {
            using T = std::decay_t< decltype( c ) >;

            if constexpr( std::is_same_v< T, EuclideanCenterDistanceConfig > )
            {
                return std::make_unique<
                    EuclideanCenterDistance< ObjectType > >(
                    c.threshold, c.weight );
            }
            else if constexpr( std::is_same_v< T, HausdorffDistanceConfig > )
            {
                return std::make_unique< HausdorffDistance< ObjectType > >(
                    c.threshold, c.weight );
            }
        },
        cfg );
}