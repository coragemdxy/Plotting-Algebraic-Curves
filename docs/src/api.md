# Public API

## Isolation and matching

- `getPointsCritical(P)` returns exact real critical projection coordinates.
- `getPointsCriticalWithIntervals(points, precise)` computes disjoint rational
  isolating intervals for those coordinates.
- `bivariateRealIsolation(P, precise=32)` isolates finite special points.
- `matchBoxesToCriticalIntervals(boxes, criticalIntervals)` matches bivariate
  boxes to exact critical coordinates.
- `criticalPointBoxes(P; boxPrecision=64, intervalPrecision=32)` runs the
  isolation and matching pipeline.

## Topology

- `getTopologySampleBoxes(P; samplesPerStrip, ...)` returns all sampled boxes.
- `topologyGraphData(P; samplesPerStrip, ...)` returns both vertices and edges.
- `connectCriticalBoxes(P; ...)` is the compatibility edge-only interface.

## Plotting and region counting

- `plotCriticalBoxEdges(edges; vertices=nothing, ...)` plots prepared graph
  data.
- `plotCriticalBoxes(P; samplesPerStrip, ...)` computes and plots a graph.
- `compactifiedTopologyGraphData(P; ...)` adds the point at infinity and all
  unbounded ends.
- `countPlaneRegions(P; ...)` returns the number of connected components of
  ``\mathbb{R}^2 \setminus \{P=0\}``.

Precision parameters are nonnegative bit counts. `boxPrecision` must be larger
than `intervalPrecision` for graph construction.
