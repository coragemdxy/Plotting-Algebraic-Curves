# Plot the graph produced by connectCriticalBoxes with Plots.jl.

using Plots

if !isdefined(@__MODULE__, :connectCriticalBoxes)
    include(joinpath(@__DIR__, "connect_critical_boxes.jl"))
end

# Return true when the current vertical direction is unsafe.  A nonconstant
# leading coefficient in y can vanish at a finite x-coordinate, allowing
# branches to escape to infinity there.  A constant leading coefficient is
# already safe even when degree(P, y) is smaller than the total degree, as for
# the parabola y - x^2.
function topologyNeedsCoordinateRotation(P)
    iszero(P) && throw(ArgumentError("P must not be the zero polynomial"))

    leadingCoefficient = leadingCoefficientInY(P)
    return degree(leadingCoefficient) > 0
end

# Put P in a safe rational coordinate system.  For a candidate integer b, use
#
#     x = u + b*v
#     y = -b*u + v.
#
# This is a rotation followed by uniform scaling.  The coefficient of v^D in
# the transformed polynomial is the highest homogeneous part of P evaluated
# at (b, 1).  A degree-D homogeneous polynomial has only finitely many bad
# directions, so the exact integer search always terminates.
function topologyRotatePolynomial(P)
    iszero(P) && throw(ArgumentError("P must not be the zero polynomial"))

    R = parent(P)
    u, v = gens(R)
    totalDegree = maximum(
        sum(exponents)
        for exponents in exponent_vectors(P)
    )

    highestPart = zero(R)
    for (coefficient, exponents) in
        zip(coefficients(P), exponent_vectors(P))
        if sum(exponents) == totalDegree
            highestPart += coefficient *
                u^exponents[1] *
                v^exponents[2]
        end
    end

    candidates = Int[0]
    for k in 1:(totalDegree + 1)
        push!(candidates, -k, k)
    end

    for candidate in candidates
        b = QQ(candidate)
        d = QQ(1)

        iszero(evaluate(highestPart, [b, d])) && continue

        originalX = d * u + b * v
        originalY = -b * u + d * v
        transformedP = evaluate(P, [originalX, originalY])

        degree(transformedP, 2) == totalDegree || continue
        leadingCoefficient = leadingCoefficientInY(transformedP)
        degree(leadingCoefficient) == 0 || continue
        iszero(leadingCoefficient) && continue

        return (
            polynomial = transformedP,
            b = b,
            d = d,
        )
    end

    error("could not find a safe rational coordinate rotation")
end

# Convert every representative endpoint and edge from the working (u, v)
# coordinates back to the original coordinates
#
#     x = d*u + b*v,
#     y = -b*u + d*v.
#
# Mapping is done before conversion to Float64 so that the exact algebraic
# x-coordinate and rational y-interval midpoint are retained as long as
# possible.
function topologyRotateEdgesBack(edges, transform)
    result = []
    b = transform.b
    d = transform.d

    for (firstBox, secondBox) in edges
        firstU = topologyBoxX(firstBox)
        firstVLeft, firstVRight = topologyBoxYInterval(firstBox)
        firstV = (firstVLeft + firstVRight) / 2
        firstPoint = (
            Float64(d * firstU + b * firstV),
            Float64(-b * firstU + d * firstV),
        )

        secondU = topologyBoxX(secondBox)
        secondVLeft, secondVRight = topologyBoxYInterval(secondBox)
        secondV = (secondVLeft + secondVRight) / 2
        secondPoint = (
            Float64(d * secondU + b * secondV),
            Float64(-b * secondU + d * secondV),
        )

        push!(result, (firstPoint, secondPoint))
    end

    return result
end


# Use the exact algebraic x-coordinate and the midpoint of the isolating
# y-interval as the displayed representative of a topology box.  Edges that
# were already mapped back to Float64 point coordinates pass through unchanged.
function criticalBoxPlotPoint(box)
    box isa Tuple{Float64, Float64} && return box

    exactX = topologyBoxX(box)
    yLeft, yRight = topologyBoxYInterval(box)
    return (Float64(exactX), Float64((yLeft + yRight) / 2))
end


#    plotCriticalBoxEdges(edges; kwargs...)

#Plot box edges returned by `connectCriticalBoxes`, or Float64 point edges that
#have already been mapped back by `topologyRotateEdgesBack`.

#Each endpoint box is represented by its exact x-coordinate (converted to
#`Float64`) and the midpoint of its rational y-isolating interval.  The result
#is a `Plots.Plot`, so callers can further customize, display, or save it.

function plotCriticalBoxEdges(
    edges;
    title::AbstractString = "Topology graph of P(x, y) = 0",
    edgeColor = :steelblue,
    vertexColor = :darkorange,
    lineWidth::Real = 2,
    markerSize::Real = 5,
)
    plt = plot(
        ;
        xlabel = "x",
        ylabel = "y",
        title = title,
        legend = false,
        aspect_ratio = :equal,
        grid = true,
    )

    vertices = Tuple{Float64, Float64}[]

    for (firstBox, secondBox) in edges
        firstPoint = criticalBoxPlotPoint(firstBox)
        secondPoint = criticalBoxPlotPoint(secondBox)

        plot!(
            plt,
            [firstPoint[1], secondPoint[1]],
            [firstPoint[2], secondPoint[2]];
            color = edgeColor,
            linewidth = lineWidth,
        )

        push!(vertices, firstPoint, secondPoint)
    end

    unique!(vertices)
    if !isempty(vertices)
        scatter!(
            plt,
            first.(vertices),
            last.(vertices);
            color = vertexColor,
            markersize = markerSize,
            markerstrokewidth = 0,
        )
    end

    return plt
end

#=
    plotCriticalBoxes(P; samplesPerStrip, output=nothing,
                      topology options..., plot options...)

Compute the edges with `connectCriticalBoxes(P)` and draw the corresponding
topology graph.  `samplesPerStrip` controls how many ordinary sample fibers
are used in every bounded and outer strip.  If the current vertical direction
has a nonconstant leading coefficient, `P` is first moved to a safe rational
coordinate system; the computed edges are mapped back before plotting.  When
`output` is a path, the figure is also saved there.
=#
function plotCriticalBoxes(
    P;
    output::Union{Nothing, AbstractString} = nothing,
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
    samplesPerStrip::Int,
    refinementStep::Int = 8,
    maxRefinements::Int = 8,
    title::AbstractString = "Topology graph of P(x, y) = 0",
    edgeColor = :steelblue,
    vertexColor = :darkorange,
    lineWidth::Real = 2,
    markerSize::Real = 5,
)
    needsRotation = topologyNeedsCoordinateRotation(P)
    workingP = P
    transform = (
        b = QQ(0),
        d = QQ(1),
    )

    if needsRotation
        rotationData = topologyRotatePolynomial(P)
        workingP = rotationData.polynomial
        transform = (
            b = rotationData.b,
            d = rotationData.d,
        )
    end

    workingEdges = connectCriticalBoxes(
        workingP;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
        samplesPerStrip = samplesPerStrip,
        refinementStep = refinementStep,
        maxRefinements = maxRefinements,
    )

    edges = topologyRotateEdgesBack(workingEdges, transform)

    plt = plotCriticalBoxEdges(
        edges;
        title = title,
        edgeColor = edgeColor,
        vertexColor = vertexColor,
        lineWidth = lineWidth,
        markerSize = markerSize,
    )


    return plt
end

function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    P = (x^2+y^2-1)^3-x^2*y^3

    output = joinpath(@__DIR__, "critical_box_graph.png")
    plotCriticalBoxes(
        P;
        samplesPerStrip = 10,
        output = output,
        boxPrecision = 32,
        intervalPrecision = 16,
    )
end
