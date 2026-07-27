# Plot the graph produced by connectCriticalBoxes with Plots.jl.

using Plots

if !isdefined(@__MODULE__, :connectCriticalBoxes)
    include(joinpath(@__DIR__, "connect_critical_boxes.jl"))
end


# Use the exact algebraic x-coordinate and the midpoint of the isolating
# y-interval as the displayed representative of a topology box.
function criticalBoxPlotPoint(box)
    exactX = topologyBoxX(box)
    yLeft, yRight = topologyBoxYInterval(box)
    return (Float64(exactX), Float64((yLeft + yRight) / 2))
end

"""
    plotCriticalBoxEdges(edges; kwargs...)

Plot the edges returned by `connectCriticalBoxes`.

Each endpoint box is represented by its exact x-coordinate (converted to
`Float64`) and the midpoint of its rational y-isolating interval.  The result
is a `Plots.Plot`, so callers can further customize, display, or save it.
"""
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

"""
    plotCriticalBoxes(P; output=nothing, topology options..., plot options...)

Compute the edges with `connectCriticalBoxes(P)` and draw the corresponding
topology graph.  When `output` is a path, the figure is also saved there.
"""
function plotCriticalBoxes(
    P;
    output::Union{Nothing, AbstractString} = nothing,
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
    refinementStep::Int = 8,
    maxRefinements::Int = 8,
    title::AbstractString = "Topology graph of P(x, y) = 0",
    edgeColor = :steelblue,
    vertexColor = :darkorange,
    lineWidth::Real = 2,
    markerSize::Real = 5,
)
    edges = connectCriticalBoxes(
        P;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
        refinementStep = refinementStep,
        maxRefinements = maxRefinements,
    )

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


# Small runnable example:
#
#   julia --project=src src/plot_critical_box_edges.jl
#
function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    P = (x^2 + y^2 - 1)

    output = joinpath(@__DIR__, "critical_box_graph.png")
    plotCriticalBoxes(
        P;
        output = output,
        boxPrecision = 32,
        intervalPrecision = 16,
    )
end
