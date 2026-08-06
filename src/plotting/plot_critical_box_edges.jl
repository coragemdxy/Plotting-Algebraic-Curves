# Plot the graph produced by connectCriticalBoxes with Plots.jl.

using Plots

if !isdefined(@__MODULE__, :connectCriticalBoxes)
    include(joinpath(@__DIR__, "..", "topology", "connect_critical_boxes.jl"))
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

# Convert one representative vertex from the working (u, v) coordinates back
# to the original coordinates
#
#     x = d*u + b*v,
#     y = -b*u + d*v.
#
# Mapping is done before conversion to Float64 so that the exact algebraic
# x-coordinate and rational y-interval midpoint are retained as long as
# possible.
function topologyRotateVertexBack(vertex, transform)
    b = transform.b
    d = transform.d

    u = topologyBoxX(vertex)
    vLeft, vRight = topologyBoxYInterval(vertex)
    v = (vLeft + vRight) / 2
    return (
        Float64(d * u + b * v),
        Float64(-b * u + d * v),
    )
end

# Convert every edge from working coordinates back to original coordinates.
function topologyRotateEdgesBack(edges, transform)
    result = []
    for (firstBox, secondBox) in edges
        push!(
            result,
            (
                topologyRotateVertexBack(firstBox, transform),
                topologyRotateVertexBack(secondBox, transform),
            ),
        )
    end
    return result
end

function topologyRotateVerticesBack(vertices, transform)
    return [
        topologyRotateVertexBack(vertex, transform)
        for vertex in vertices
    ]
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


"""
    plotCriticalBoxEdges(edges; vertices=nothing, kwargs...)

Plot box edges returned by `connectCriticalBoxes`, or Float64 point edges that
have already been mapped back by `topologyRotateEdgesBack`.  Passing the
complete `vertices` collection also plots degree-zero vertices, including
isolated real points.

Each box is represented by its exact x-coordinate (converted to `Float64`) and
the midpoint of its rational y-isolating interval.  The result is a
`Plots.Plot`, so callers can further customize, display, or save it.
"""

function plotCriticalBoxEdges(
    edges;
    vertices = nothing,
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

    plottedVertices = Tuple{Float64, Float64}[]

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

        push!(plottedVertices, firstPoint, secondPoint)
    end

    if vertices !== nothing
        for vertex in vertices
            push!(plottedVertices, criticalBoxPlotPoint(vertex))
        end
    end

    unique!(plottedVertices)
    if !isempty(plottedVertices)
        scatter!(
            plt,
            first.(plottedVertices),
            last.(plottedVertices);
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

Compute the vertices and edges with `topologyGraphData(P)` and draw the
corresponding topology graph.  Retaining the explicit vertex collection makes
isolated real points visible even though they have no incident edge.
`samplesPerStrip` controls how many ordinary sample fibers are used in every
bounded and outer strip.  If the current vertical direction has a nonconstant
leading coefficient, `P` is first moved to a safe rational coordinate system;
the computed vertices and edges are mapped back before plotting.  When
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

    workingGraph = topologyGraphData(
        workingP;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
        samplesPerStrip = samplesPerStrip,
        refinementStep = refinementStep,
        maxRefinements = maxRefinements,
    )

    edges = topologyRotateEdgesBack(workingGraph.edges, transform)
    vertices =
        topologyRotateVerticesBack(workingGraph.vertices, transform)

    plt = plotCriticalBoxEdges(
        edges;
        vertices = vertices,
        title = title,
        edgeColor = edgeColor,
        vertexColor = vertexColor,
        lineWidth = lineWidth,
        markerSize = markerSize,
    )

    if output !== nothing
        savefig(plt, output)
    end
    return plt
end
