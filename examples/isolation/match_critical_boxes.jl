using PlotCurveTopologyForBivariatePolynomial
using Nemo

function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 1) * (x^2 + y^2 - 4) * x
    box = criticalPointBoxes(f)
    return box
end
