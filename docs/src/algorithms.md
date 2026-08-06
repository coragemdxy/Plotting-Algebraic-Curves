# Algorithms

The active implementation follows this pipeline:

1. Remove repeated factors and split vertical components from the primitive
   nonvertical part.
2. Compute critical projection coordinates by exact elimination.
3. Isolate special points in rational boxes and match them to exact algebraic
   critical coordinates.
4. Add ordinary sample fibers and certify branch incidence near every critical
   fiber with horizontal separating boundaries.
5. Build the sampled topology graph, retaining isolated vertices explicitly.
6. For region counting, choose a safe projection direction when necessary,
   compactify the graph with one vertex at infinity, and apply Euler's formula.

The historical rational-only implementation is retained under
`prototype/rational_method/`. Supporting mathematical proofs are stored in
`docs/proofs/`, and the internship report is stored in `docs/report/`.
