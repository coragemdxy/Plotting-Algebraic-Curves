# PlotCurveTopologyForBivariatePolynomial.jl

[![CI](https://github.com/coragemdxy/PlotCurveTopologyForBivariatePolynomial.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/coragemdxy/PlotCurveTopologyForBivariatePolynomial.jl/actions/workflows/CI.yml)

`PlotCurveTopologyForBivariatePolynomial.jl` computes certified sampled topology graphs for real
plane algebraic curves defined over `QQ`. It can isolate critical points,
connect branches across critical fibers, plot a representative graph, and count
the connected regions of the curve complement.

The project is an exact computer-algebra research implementation. Topological
decisions use rational numbers, exact algebraic numbers, and certified interval
operations; floating-point values are used only to display representative
points.

## Installation

Until the package is registered in Julia's General registry, install it from
GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/coragemdxy/PlotCurveTopologyForBivariatePolynomial.jl.git")
```

The supported Julia version is 1.10 or later. The direct dependencies are
AlgebraicSolving, Nemo, and Plots.

## Quick start

```julia
using PlotCurveTopologyForBivariatePolynomial
using Nemo

R, (x, y) = polynomial_ring(QQ, ["x", "y"])
P = (x^2 + y^2 - 1) * (x^2 + y^2 - 4)

graph = topologyGraphData(
    P;
    samplesPerStrip = 1,
    boxPrecision = 64,
    intervalPrecision = 32,
)

graph.vertices
graph.edges

countPlaneRegions(P) # 3
```

To draw the sampled topology graph:

```julia
plotCriticalBoxes(
    P;
    samplesPerStrip = 4,
    boxPrecision = 64,
    intervalPrecision = 32,
)
```

## Main interfaces

- `bivariateRealIsolation` isolates finite special points in rational boxes.
- `criticalPointBoxes` associates those boxes with exact critical coordinates.
- `topologyGraphData` returns the complete sampled vertex and edge collections.
- `connectCriticalBoxes` retains the historical edge-only interface.
- `plotCriticalBoxes` plots the graph in the original coordinates.
- `compactifiedTopologyGraphData` adds the point at infinity and unbounded ends.
- `countPlaneRegions` counts components of `R² \ {P = 0}`.

See [the API guide](docs/src/api.md) for parameters and supporting interfaces.

## Repository layout

- `src/` contains the registered-package implementation.
- `test/` contains exact regression tests for isolation, topology, and region
  counting.
- `examples/` contains the original interactive `main` examples, moved out of
  package source without changing their implementations.
- `prototype/rational_method/` preserves the earlier rational-only Method 1.
- `docs/proofs/` contains mathematical LaTeX demonstrations.
- `docs/report/` contains the internship report.

## Development

Instantiate the package environment and run the tests:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

See [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

## Scope and limitations

- Input is expected to be a bivariate polynomial over `QQ`.
- Exact elimination and algebraic-number matching may become expensive as
  degree and coefficient size grow.
- Plot segments are visual representatives of certified graph incidence; they
  are not linear approximations used for topology decisions.
- The rational Method 1 is retained for historical and pedagogical purposes but
  is not loaded by the package.

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff). The full
mathematical and implementation background is available in the
[internship report](docs/report/internship_report.pdf).

## License

This project is available under the [MIT License](LICENSE).
