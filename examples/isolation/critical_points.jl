using PlotCurveTopologyForBivariatePolynomial
using Nemo

function main()
    R,(x,y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 3)
    println(f)
    s = getPointsCritical(f)
    s1 = getPointsCriticalWithIntervals(s,32)
    println(s1)
end
