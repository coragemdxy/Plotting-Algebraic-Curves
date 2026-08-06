using PlotCurveTopologyForBivariatePolynomial
using Nemo

function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 1) * (x^2 + y^2 - 4) * x
    boxes = bivariateRealIsolation(f, 16)
    println("f = ", f)
    display(boxes)
end
