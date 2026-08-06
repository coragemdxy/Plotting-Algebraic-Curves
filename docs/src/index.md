# PlotCurveTopologyForBivariatePolynomial.jl

`PlotCurveTopologyForBivariatePolynomial.jl` computes a certified sampled topology graph for real
plane algebraic curves defined by polynomials in ``\mathbb{Q}[x,y]``. The same
graph can be plotted or compactified to count the connected regions of the
curve complement.

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

countPlaneRegions(P) # 3
```

The topology decisions use exact rational and algebraic arithmetic. Floating
point conversion is confined to plotting representative points.
