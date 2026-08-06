using PlotCurveTopologyForBivariatePolynomial
using Nemo

function main(samplesPerStrip::Int)
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    P = (x^2 + y^2 - 1) * x

    edges = connectCriticalBoxes(
        P;
        samplesPerStrip = samplesPerStrip,
        boxPrecision = 32,
        intervalPrecision = 16,
    )

    println(edges)
    println(length(edges))
end
