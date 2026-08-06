#=
Count the connected regions of R^2 \ {P = 0}.

The curve graph constructed by connect_critical_boxes.jl contains all finite
critical and ordinary-fiber vertices, but its unbounded branches stop at the
outer sample fibers.  This file compactifies the plane to a sphere: it adds
one vertex at infinity and joins every unbounded branch end to that vertex.

For a finite graph embedded in the sphere, with V vertices, E edges and C
connected components, Euler's formula gives

    number of regions = E - V + C + 1.

The main interface is

    countPlaneRegions(P; samplesPerStrip = 1, ...)

and compactifiedTopologyGraphData(P; ...) exposes the compactified graph used
by the count.
=#

using Nemo

if !isdefined(@__MODULE__, :topologyGraphData)
    include(joinpath(@__DIR__, "connect_critical_boxes.jl"))
end


# The extra graph vertex representing the point at infinity of the sphere.
# Keep the guard so that the file can be included repeatedly in one session.
if !isdefined(@__MODULE__, :PlaneInfinityVertex)
    struct PlaneInfinityVertex end
end


#=
Return true when the current vertical projection can have branches escaping
to infinity above a finite x-coordinate.  A constant leading coefficient in
y rules this out.
=#
function planeRegionNeedsCoordinateRotation(P)
    iszero(P) && throw(ArgumentError("P must not be the zero polynomial"))
    return degree(leadingCoefficientInY(P)) > 0
end


#=
Put P in coordinates in which its leading coefficient in the second
variable is a nonzero constant.  For an integer b, use

     x = u + b*v
     y = -b*u + v.

The determinant is 1 + b^2, so this is an invertible linear map of R^2 and
does not change the number of connected regions in the complement.
=#
function planeRegionRotatePolynomial(P)
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

    # The polynomial highestPart(b, 1) has degree at most totalDegree in b,
    # hence only finitely many bad integer directions can occur.
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
        leadingCoefficient = leadingCoefficientInY(transformedP)

        degree(transformedP, 2) == totalDegree || continue
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


# Return P in a safe coordinate system and retain the linear transformation.
function planeRegionSafePolynomial(P)
    Psf = bivariateSquarefreePart(P)

    if planeRegionNeedsCoordinateRotation(Psf)
        return planeRegionRotatePolynomial(Psf)
    end

    return (
        polynomial = Psf,
        b = QQ(0),
        d = QQ(1),
    )
end


#= 
The first fiber in the left outer strip and the last fiber in the right
outer strip are the outermost sampled fibers.  In safe coordinates, every
root on one of these fibers represents exactly one unbounded curve end.
=#
function planeRegionUnboundedEndBoxes(data)
    isempty(data.ordinaryBoxGroupsByStrip) && return []

    leftOuterStrip = first(data.ordinaryBoxGroupsByStrip)
    rightOuterStrip = last(data.ordinaryBoxGroupsByStrip)
    isempty(leftOuterStrip) && return []
    isempty(rightOuterStrip) && return []

    leftEndBoxes = first(leftOuterStrip)
    rightEndBoxes = last(rightOuterStrip)
    return vcat(leftEndBoxes, rightEndBoxes)
end


#=
Construct the graph of the compactified zero set in S^2.  The point at
infinity is always retained.  For a compact curve it is an isolated graph
vertex; for an unbounded curve it is incident to all unbounded branch ends.
=#
function compactifiedTopologyGraphData(
    P;
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
    samplesPerStrip::Int = 1,
    refinementStep::Int = 8,
    maxRefinements::Int = 8,
)
    requireBivariateQQPolynomial(P)
    iszero(P) && throw(ArgumentError(
        "the zero polynomial does not define a finite curve graph",
    ))
    boxPrecision >= 0 ||
        throw(ArgumentError("boxPrecision must be nonnegative"))
    intervalPrecision >= 0 ||
        throw(ArgumentError("intervalPrecision must be nonnegative"))
    refinementStep > 0 ||
        throw(ArgumentError("refinementStep must be positive"))
    maxRefinements >= 0 ||
        throw(ArgumentError("maxRefinements must be nonnegative"))
    samplesPerStrip >= 1 ||
        throw(ArgumentError("samplesPerStrip must be at least 1"))
    boxPrecision > intervalPrecision ||
        throw(ArgumentError(
            "boxPrecision must be larger than intervalPrecision",
        ))

    transform = planeRegionSafePolynomial(P)
    workingP = transform.polynomial

    data = topologyPrepareBoxData(
        workingP;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
        samplesPerStrip = samplesPerStrip,
    )
    finiteEdges = topologyEdgesFromPreparedData(
        data;
        intervalPrecision = intervalPrecision,
        refinementStep = refinementStep,
        maxRefinements = maxRefinements,
    )

    vertices = Any[data.allBoxes...]
    edges = Any[finiteEdges...]
    infinityVertex = PlaneInfinityVertex()
    push!(vertices, infinityVertex)

    unboundedEndBoxes = planeRegionUnboundedEndBoxes(data)
    for box in unboundedEndBoxes
        push!(edges, (box, infinityVertex))
    end

    return (
        vertices = vertices,
        edges = edges,
        infinityVertex = infinityVertex,
        unboundedEndBoxes = unboundedEndBoxes,
        workingPolynomial = workingP,
        transform = (
            b = transform.b,
            d = transform.d,
        ),
    )
end


# Find a representative of one disjoint-set class and compress its path.
function planeRegionFindRoot!(parents, vertex::Int)
    current = vertex
    while parents[current] != current
        current = parents[current]
    end

    root = current
    current = vertex
    while parents[current] != current
        next = parents[current]
        parents[current] = root
        current = next
    end
    return root
end


# Merge two disjoint-set classes by rank.
function planeRegionUnion!(parents, ranks, firstVertex::Int, secondVertex::Int)
    firstRoot = planeRegionFindRoot!(parents, firstVertex)
    secondRoot = planeRegionFindRoot!(parents, secondVertex)
    firstRoot == secondRoot && return

    if ranks[firstRoot] < ranks[secondRoot]
        parents[firstRoot] = secondRoot
    elseif ranks[secondRoot] < ranks[firstRoot]
        parents[secondRoot] = firstRoot
    else
        parents[secondRoot] = firstRoot
        ranks[firstRoot] += 1
    end
end


#= 
Count connected components of a finite multigraph.  Edges are deliberately
not deduplicated: two edges with the same endpoints can represent distinct
curve arcs and both must contribute to Euler's formula.
=#
function planeRegionGraphComponentCount(vertices, edges)
    numberOfVertices = length(vertices)
    numberOfVertices == 0 && return 0

    vertexIndices = Dict{Any, Int}()
    for (i, vertex) in enumerate(vertices)
        haskey(vertexIndices, vertex) &&
            throw(ArgumentError("the graph contains a duplicate vertex"))
        vertexIndices[vertex] = i
    end

    parents = collect(1:numberOfVertices)
    ranks = zeros(Int, numberOfVertices)

    for (firstVertex, secondVertex) in edges
        haskey(vertexIndices, firstVertex) ||
            throw(ArgumentError("an edge endpoint is missing from vertices"))
        haskey(vertexIndices, secondVertex) ||
            throw(ArgumentError("an edge endpoint is missing from vertices"))

        planeRegionUnion!(
            parents,
            ranks,
            vertexIndices[firstVertex],
            vertexIndices[secondVertex],
        )
    end

    components = Set{Int}()
    for i in 1:numberOfVertices
        push!(components, planeRegionFindRoot!(parents, i))
    end
    return length(components)
end


#=
    countPlaneRegions(P; samplesPerStrip=1,
                      boxPrecision=64, intervalPrecision=32,
                      refinementStep=8, maxRefinements=8)

Return the number of connected components of `R^2 \\ {P = 0}`.

The polynomial must belong to `QQ[x,y]`.  Repeated factors are removed, a
safe rational coordinate system is chosen when necessary, and the topology
graph is compactified by adding its point at infinity.  The final count uses
`E - V + C + 1` on the sphere.
=#
function countPlaneRegions(
    P;
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
    samplesPerStrip::Int = 1,
    refinementStep::Int = 8,
    maxRefinements::Int = 8,
)
    # If P is identically zero, its zero set is the entire plane and the
    # complement has no connected component.
    iszero(P) && return 0
    requireBivariateQQPolynomial(P)

    graph = compactifiedTopologyGraphData(
        P;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
        samplesPerStrip = samplesPerStrip,
        refinementStep = refinementStep,
        maxRefinements = maxRefinements,
    )

    numberOfVertices = length(graph.vertices)
    numberOfEdges = length(graph.edges)
    numberOfComponents = planeRegionGraphComponentCount(
        graph.vertices,
        graph.edges,
    )

    return numberOfEdges - numberOfVertices + numberOfComponents + 1
end


# A small example for interactive use.  This function is not run when the file is included; call main() explicitly.
