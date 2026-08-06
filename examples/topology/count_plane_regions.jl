using PlotCurveTopologyForBivariatePolynomial
using Nemo

function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    P = x*y-1

    numberOfRegions = countPlaneRegions(
        P;
        samplesPerStrip = 1,
        boxPrecision = 32,
        intervalPrecision = 16,
    )

    println("P = ", P)
    println("number of connected regions = ", numberOfRegions)
    return numberOfRegions
end
