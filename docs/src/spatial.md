```@meta
CurrentModule = BranchingProcesses
```

# Spatial Mapping

Functions for adding 2D or 3D spatial locations to branching process trees to model cell proliferation where daughter cells stay in close physical proximity after division.

## Spatial Tree Type

```@docs
SpatialTree
```

## Position Assignment

```@docs
assign_spatial_positions
```

## Position Refinement

```@docs
relax_positions!
normalize_positions!
```

## Position and Value Extraction

```@docs
get_leaf_positions
get_positions_at_time
get_values_at_positions
```

## Convenience Functions

```@docs
create_spatial_layout
```

## Visualization Helpers

```@docs
spatial_heatmap_data
get_time_series_frames
```
