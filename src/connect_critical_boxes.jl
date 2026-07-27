#=
Connect critical-fiber boxes and ordinary sample-fiber boxes.

The main function is

    connectCriticalBoxes(P)

It returns an array of edges.  Every edge has the form

    (firstBox, secondBox)

Every box keeps the format used by match_critical_boxes.jl:

    ((exactCriticalX, criticalXInterval), yInterval)

Besides the true critical x-coordinates, the program adds one rational sample
x-coordinate in every open strip and one sample coordinate on each outer
side.  Consequently an unbounded branch is represented up to an outer sample
box instead of disappearing from the edge list.
=#

using Nemo

# Loading match_critical_boxes.jl also loads Qbar_critical_points.jl and
# bivariate_real_isolation.jl.  The guard makes this file safe to include after
# those files have already been loaded.
if !isdefined(@__MODULE__, :criticalPointBoxes)
    include(joinpath(@__DIR__, "match_critical_boxes.jl"))
end



# Exact critical x-coordinate of a matched box.
topologyBoxX(box) = box[1][1]

# Rational y-interval of a matched box.
topologyBoxYInterval(box) = box[2]

# A rational number used only for sorting boxes on the same vertical line.
function topologyBoxYMiddle(box)
    yLeft, yRight = topologyBoxYInterval(box)
    return (yLeft + yRight) / 2
end

# Return all boxes on x = alpha, ordered from bottom to top.
function topologyBoxesAtX(boxes, alpha)
    result = []

    for box in boxes
        if topologyBoxX(box) == alpha
            push!(result, box)
        end
    end

    sort!(result, by = topologyBoxYMiddle)
    return result
end

# Return the exact real roots of H(beta, y), ordered from bottom to top.
# beta is rational, so this is an ordinary fiber and no QQBar-to-QQ conversion
# is needed.
function topologyRootsOnVerticalFiber(H, beta)
    R = parent(H)
    _, y = gens(R)

    polynomialInY = evaluate(H, [QQ(beta), y])
    univariatePolynomial = toUnivariate(polynomialInY, 2)

    univariatePolynomial === nothing &&
        error("internal error: H(beta, y) still depends on x")

    # This should not happen after verticalDecomposition, because H is
    # primitive with respect to y.  Keeping the check gives a clearer error if
    # the function is called with unsuitable input.
    iszero(univariatePolynomial) &&
        throw(ArgumentError("H(beta, y) is the zero polynomial at beta = $beta"))

    rootsAtBeta = realRootsAsAlgebraic(univariatePolynomial)
    sort!(rootsAtBeta)
    return rootsAtBeta
end

# Choose the ordinary rational x-coordinates used as extra sample fibers.
# For critical intervals I_1, ..., I_s the result is
#     beta_0, beta_1, ..., beta_s,
# where beta_0 is left of I_1, beta_s is right of I_s, and beta_i lies
# strictly between I_i and I_(i+1).
function topologyOrdinarySampleXs(criticalIntervals)
    if isempty(criticalIntervals)
        # A curve such as y = 0 has no critical x-coordinate.  Two ordinary
        # fibers are still needed in order to produce one representative edge.
        return [QQ(-1), QQ(1)]
    end

    sampleXs = []

    firstInterval = criticalIntervals[1][2]
    push!(sampleXs, firstInterval[1] - 1)

    for i in 1:(length(criticalIntervals) - 1)
        currentRight = criticalIntervals[i][2][2]
        nextLeft = criticalIntervals[i + 1][2][1]

        currentRight < nextLeft ||
            error("critical x-intervals must be disjoint")

        push!(sampleXs, (currentRight + nextLeft) / 2)
    end

    lastInterval = criticalIntervals[end][2]
    push!(sampleXs, lastInterval[2] + 1)

    return sampleXs
end

# Isolate all points of H = 0 on the ordinary rational fiber x = beta and use
# matchBoxesToCriticalIntervals to give them exactly the same representation
# as the true critical boxes.
function topologyBoxesOnOrdinaryFiber(H, beta, boxPrecision::Int)
    R = parent(H)
    x, _ = gens(R)

    points = zeroDimensionalPoints(H, x - QQ(beta), boxPrecision)
    rawBoxes = [point.box for point in points]

    isempty(rawBoxes) && return []

    Qb = algebraic_closure(QQ)
    exactBeta = Qb(beta)

    # isolateAlgebraicNumber sometimes returns a very small one-sided interval
    # even when beta is rational.  Use the union of those x-intervals as the
    # matching interval instead of assuming that every interval is (beta,beta).
    sampleLeft = minimum(box[1][1] for box in rawBoxes)
    sampleRight = maximum(box[1][2] for box in rawBoxes)
    sampleInterval = (sampleLeft, sampleRight)

    return matchBoxesToCriticalIntervals(
        rawBoxes,
        [(exactBeta, sampleInterval)],
    )
end

# Sort arbitrary sample boxes first by exact x and then by y.
function topologyBoxIsLess(firstBox, secondBox)
    firstX = topologyBoxX(firstBox)
    secondX = topologyBoxX(secondBox)
    firstX == secondX || return firstX < secondX

    return topologyBoxYMiddle(firstBox) < topologyBoxYMiddle(secondBox)
end


# Substitute y = height and return the resulting polynomial in x.
function topologyHorizontalSlice(H, height)
    R = parent(H)
    x, _ = gens(R)

    polynomialInX = evaluate(H, [x, QQ(height)])
    univariatePolynomial = toUnivariate(polynomialInX, 1)

    univariatePolynomial === nothing &&
        error("internal error: H(x, height) still depends on y")

    return univariatePolynomial
end

# Choose a rational horizontal boundary strictly between low and high.
# A boundary which is itself a horizontal component of H is rejected.
function topologyChooseHorizontalBoundary(H, low, high)
    low < high ||
        throw(ArgumentError("a horizontal boundary needs low < high"))

    # H has at most degree(H, y) distinct horizontal components.  Trying a few
    # more candidates guarantees that at least one candidate is not such a
    # component.
    numberOfTries = Int(degree(H, 2)) + 3

    for k in 1:numberOfTries
        fraction = QQ(k, numberOfTries + 1)
        candidate = low + fraction * (high - low)

        if !iszero(topologyHorizontalSlice(H, candidate))
            return candidate
        end
    end

    error("could not choose a horizontal boundary between $low and $high")
end

# Build disjoint open y-windows around the boxes on one critical fiber.
# The returned boundaries have the form
#
#     b_1 < b_2 < ... < b_(n+1),
#
# and box j lies in the window (b_j, b_(j+1)).
function topologyBuildYWindows(H, boxesAtAlpha)
    isempty(boxesAtAlpha) && return ([], [])

    orderedBoxes = sort(collect(boxesAtAlpha), by = topologyBoxYMiddle)
    numberOfBoxes = length(orderedBoxes)

    # Check that the isolating y-intervals are pairwise disjoint.
    for j in 1:(numberOfBoxes - 1)
        currentRight = topologyBoxYInterval(orderedBoxes[j])[2]
        nextLeft = topologyBoxYInterval(orderedBoxes[j + 1])[1]

        currentRight < nextLeft ||
            throw(ArgumentError(
                "the y-intervals of two boxes on the same critical fiber " *
                "overlap; increase boxPrecision",
            ))
    end

    boundaries = []

    firstLeft = topologyBoxYInterval(first(orderedBoxes))[1]
    push!(
        boundaries,
        topologyChooseHorizontalBoundary(H, firstLeft - 1, firstLeft),
    )

    for j in 1:(numberOfBoxes - 1)
        currentRight = topologyBoxYInterval(orderedBoxes[j])[2]
        nextLeft = topologyBoxYInterval(orderedBoxes[j + 1])[1]
        push!(
            boundaries,
            topologyChooseHorizontalBoundary(H, currentRight, nextLeft),
        )
    end

    lastRight = topologyBoxYInterval(last(orderedBoxes))[2]
    push!(
        boundaries,
        topologyChooseHorizontalBoundary(H, lastRight, lastRight + 1),
    )

    windows = []
    for j in 1:numberOfBoxes
        push!(windows, (boundaries[j], boundaries[j + 1]))
    end

    return windows, boundaries
end

# Compute the exact x-coordinates where the curve meets each horizontal
# boundary.  These roots can be reused while the critical x-interval is being
# refined.
function topologyRootsOnHorizontalBoundaries(H, boundaries)
    result = []

    for boundary in boundaries
        polynomialInX = topologyHorizontalSlice(H, boundary)
        rootsOnBoundary = realRootsAsAlgebraic(polynomialInX)
        sort!(rootsOnBoundary)
        push!(result, rootsOnBoundary)
    end

    return result
end

# True when no horizontal boundary meets the curve above the x-interval
# [xLeft, xRight].  This is the certification step which prevents a branch
# from leaving one y-window and entering another one near the critical fiber.
function topologyBoundariesAreClear(boundaryRoots, xLeft, xRight)
    for rootsOnOneBoundary in boundaryRoots
        for root in rootsOnOneBoundary
            Qb = parent(root)
            if Qb(xLeft) <= root <= Qb(xRight)
                return false
            end
        end
    end

    return true
end

# For each ordered root, return the index of its unique y-window.  Zero means
# that the root has no finite box endpoint.  This can happen at a root of the
# leading coefficient, where a branch tends to infinity.
function topologyAssignRootsToWindows(rootsAtX, windows)
    assignments = zeros(Int, length(rootsAtX))

    for k in eachindex(rootsAtX)
        root = rootsAtX[k]
        Qb = parent(root)

        for j in eachindex(windows)
            lower, upper = windows[j]
            if Qb(lower) < root < Qb(upper)
                assignments[k] = j
                break
            end
        end
    end

    return assignments
end

# Find the isolating interval belonging to one exact critical x-coordinate.
function topologyFindCriticalInterval(criticalIntervals, alpha)
    for (point, interval) in criticalIntervals
        if point == alpha
            return interval
        end
    end

    error("could not find an isolating interval for critical point $alpha")
end

# Determine which critical box each branch approaches from the left and from
# the right.  The critical interval is refined until all horizontal boundaries
# are certified to be free of curve points in the local slab.
function topologyIncidenceAtCriticalX(
    H,
    criticalPoints,
    alpha,
    boxesAtAlpha;
    intervalPrecision::Int,
    refinementStep::Int,
    maxRefinements::Int,
)
    orderedBoxes = sort(collect(boxesAtAlpha), by = topologyBoxYMiddle)

    # A critical value caused only by the leading coefficient may have no
    # finite point above it.  Its neighboring branches are unbounded, so none
    # of them is assigned to a finite box.
    if isempty(orderedBoxes)
        intervalData =
            getPointsCriticalWithIntervals(criticalPoints, intervalPrecision)
        xLeft, xRight = topologyFindCriticalInterval(intervalData, alpha)
        leftRoots = topologyRootsOnVerticalFiber(H, xLeft)
        rightRoots = topologyRootsOnVerticalFiber(H, xRight)

        return (
            boxes = orderedBoxes,
            leftRoots = leftRoots,
            rightRoots = rightRoots,
            leftAssignments = zeros(Int, length(leftRoots)),
            rightAssignments = zeros(Int, length(rightRoots)),
        )
    end

    windows, boundaries = topologyBuildYWindows(H, orderedBoxes)
    boundaryRoots = topologyRootsOnHorizontalBoundaries(H, boundaries)

    for attempt in 0:maxRefinements
        currentPrecision = intervalPrecision + attempt * refinementStep
        intervalData =
            getPointsCriticalWithIntervals(criticalPoints, currentPrecision)
        xLeft, xRight = topologyFindCriticalInterval(intervalData, alpha)

        if topologyBoundariesAreClear(boundaryRoots, xLeft, xRight)
            leftRoots = topologyRootsOnVerticalFiber(H, xLeft)
            rightRoots = topologyRootsOnVerticalFiber(H, xRight)

            return (
                boxes = orderedBoxes,
                leftRoots = leftRoots,
                rightRoots = rightRoots,
                leftAssignments =
                    topologyAssignRootsToWindows(leftRoots, windows),
                rightAssignments =
                    topologyAssignRootsToWindows(rightRoots, windows),
            )
        end
    end

    throw(ArgumentError(
        "could not certify the local boxes above x = $alpha; " *
        "increase maxRefinements or boxPrecision",
    ))
end

# Put the box with smaller x first.  On a vertical component, put the lower
# box first.
function topologyOrderOneEdge(firstBox, secondBox)
    firstX = topologyBoxX(firstBox)
    secondX = topologyBoxX(secondBox)

    if firstX < secondX
        return (firstBox, secondBox)
    elseif secondX < firstX
        return (secondBox, firstBox)
    elseif topologyBoxYMiddle(firstBox) <= topologyBoxYMiddle(secondBox)
        return (firstBox, secondBox)
    else
        return (secondBox, firstBox)
    end
end

# Comparison used to sort the final edge list by the first critical point.
function topologyEdgeIsLess(firstEdge, secondEdge)
    firstBox1, firstBox2 = firstEdge
    secondBox1, secondBox2 = secondEdge

    x11 = topologyBoxX(firstBox1)
    x21 = topologyBoxX(secondBox1)
    x11 == x21 || return x11 < x21

    y11 = topologyBoxYMiddle(firstBox1)
    y21 = topologyBoxYMiddle(secondBox1)
    y11 == y21 || return y11 < y21

    x12 = topologyBoxX(firstBox2)
    x22 = topologyBoxX(secondBox2)
    x12 == x22 || return x12 < x22

    return topologyBoxYMiddle(firstBox2) < topologyBoxYMiddle(secondBox2)
end

# Add edges along every vertical component V(x) = 0.  Consecutive intersection
# boxes are connected from bottom to top.  The two unbounded rays are not box-
# to-box edges, so they are not included in the requested return format.
function topologyAddVerticalEdges!(edges, V, criticalPoints, boxGroups)
    for i in eachindex(criticalPoints)
        alpha = criticalPoints[i]

        if iszero(V(alpha))
            boxesOnComponent = boxGroups[i]

            for j in 1:(length(boxesOnComponent) - 1)
                push!(
                    edges,
                    topologyOrderOneEdge(
                        boxesOnComponent[j],
                        boxesOnComponent[j + 1],
                    ),
                )
            end
        end
    end
end

# Prepare all fibers and boxes used by the graph construction.  This helper is
# kept separate so that getTopologySampleBoxes and connectCriticalBoxes use
# exactly the same sample points.
function topologyPrepareBoxData(
    P;
    boxPrecision::Int,
    intervalPrecision::Int,
)
    Psf = bivariateSquarefreePart(P)
    V, H = verticalDecomposition(Psf)

    # A constant H means that P consists only of vertical components.  The
    # present box format describes isolated points on vertical fibers, not a
    # whole vertical line, so there are no finite sample boxes to prepare.
    if degree(H, 2) <= 0
        return (
            Psf = Psf,
            V = V,
            H = H,
            criticalPoints = [],
            criticalIntervals = [],
            criticalBoxes = [],
            criticalBoxGroups = [],
            ordinarySampleXs = [QQ(-1), QQ(1)],
            ordinaryBoxGroups = [[], []],
            allBoxes = [],
        )
    end

    criticalPoints = getPointsCritical(Psf)
    criticalIntervals = getPointsCriticalWithIntervals(
        criticalPoints,
        intervalPrecision,
    )

    criticalBoxes = if isempty(criticalPoints)
        []
    else
        criticalPointBoxes(
            Psf;
            boxPrecision = boxPrecision,
            intervalPrecision = intervalPrecision,
        )
    end

    criticalBoxGroups = []
    for alpha in criticalPoints
        push!(criticalBoxGroups, topologyBoxesAtX(criticalBoxes, alpha))
    end

    ordinarySampleXs = topologyOrdinarySampleXs(criticalIntervals)
    ordinaryBoxGroups = []

    for beta in ordinarySampleXs
        push!(
            ordinaryBoxGroups,
            topologyBoxesOnOrdinaryFiber(H, beta, boxPrecision),
        )
    end

    allBoxes = []
    append!(allBoxes, criticalBoxes)
    for group in ordinaryBoxGroups
        append!(allBoxes, group)
    end
    sort!(allBoxes, lt = topologyBoxIsLess)

    return (
        Psf = Psf,
        V = V,
        H = H,
        criticalPoints = criticalPoints,
        criticalIntervals = criticalIntervals,
        criticalBoxes = criticalBoxes,
        criticalBoxGroups = criticalBoxGroups,
        ordinarySampleXs = ordinarySampleXs,
        ordinaryBoxGroups = ordinaryBoxGroups,
        allBoxes = allBoxes,
    )
end

"""
    getTopologySampleBoxes(P; boxPrecision=64, intervalPrecision=32)

Return all boxes used by `connectCriticalBoxes`: true critical-fiber boxes,
one ordinary sample fiber in every strip, and the two outer sample fibers.
The boxes are sorted by x and then by y.
"""
function getTopologySampleBoxes(
    P;
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
)
    boxPrecision > intervalPrecision ||
        throw(ArgumentError(
            "boxPrecision must be larger than intervalPrecision",
        ))

    data = topologyPrepareBoxData(
        P;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
    )
    return data.allBoxes
end

"""
    connectCriticalBoxes(P; boxPrecision=64, intervalPrecision=32,
                         refinementStep=8, maxRefinements=8)

Return the box-to-box edges of the sampled real curve `P = 0`.

Each result is `(firstBox, secondBox)`.  The first box has smaller critical
x-coordinate, and the full array is sorted by that first critical point.
The vertices include ordinary midpoint and outer sample boxes, so branches
leading away from the leftmost or rightmost critical fiber are retained up to
an outer sample box.
"""
function connectCriticalBoxes(
    P;
    boxPrecision::Int = 64,
    intervalPrecision::Int = 32,
    refinementStep::Int = 8,
    maxRefinements::Int = 8,
)
    boxPrecision >= 0 ||
        throw(ArgumentError("boxPrecision must be nonnegative"))
    intervalPrecision >= 0 ||
        throw(ArgumentError("intervalPrecision must be nonnegative"))
    refinementStep > 0 ||
        throw(ArgumentError("refinementStep must be positive"))
    maxRefinements >= 0 ||
        throw(ArgumentError("maxRefinements must be nonnegative"))
    boxPrecision > intervalPrecision ||
        throw(ArgumentError(
            "boxPrecision must be larger than intervalPrecision so that " *
            "critical boxes fit inside their matched x-intervals",
        ))

    data = topologyPrepareBoxData(
        P;
        boxPrecision = boxPrecision,
        intervalPrecision = intervalPrecision,
    )

    V = data.V
    H = data.H
    criticalPoints = data.criticalPoints
    boxGroups = data.criticalBoxGroups
    ordinaryBoxGroups = data.ordinaryBoxGroups

    degree(H, 2) <= 0 && return []

    edges = []

    # If there is no critical x-coordinate, root order is constant on the
    # whole real line.  Connect the two outer ordinary fibers directly.
    if isempty(criticalPoints)
        leftBoxes = ordinaryBoxGroups[1]
        rightBoxes = ordinaryBoxGroups[2]

        length(leftBoxes) == length(rightBoxes) ||
            error("the number of roots changed although there is no critical x")

        for k in eachindex(leftBoxes)
            push!(edges, topologyOrderOneEdge(leftBoxes[k], rightBoxes[k]))
        end

        sort!(edges, lt = topologyEdgeIsLess)
        return edges
    end

    # Compute the certified left/right branch incidences of every critical
    # fiber.
    incidences = []
    for i in eachindex(criticalPoints)
        push!(
            incidences,
            topologyIncidenceAtCriticalX(
                H,
                criticalPoints,
                criticalPoints[i],
                boxGroups[i];
                intervalPrecision = intervalPrecision,
                refinementStep = refinementStep,
                maxRefinements = maxRefinements,
            ),
        )
    end

    # Every critical fiber has one ordinary sample fiber on its left and one
    # on its right.  Root order is constant between the ordinary fiber and the
    # critical fiber, so root k on both nearby ordinary lines is the same
    # branch.
    for i in eachindex(criticalPoints)
        leftOrdinaryBoxes = ordinaryBoxGroups[i]
        rightOrdinaryBoxes = ordinaryBoxGroups[i + 1]

        length(leftOrdinaryBoxes) == length(incidences[i].leftRoots) ||
            error(
                "the number of roots is inconsistent on the left of " *
                "$(criticalPoints[i])",
            )

        length(rightOrdinaryBoxes) == length(incidences[i].rightRoots) ||
            error(
                "the number of roots is inconsistent on the right of " *
                "$(criticalPoints[i])",
            )

        for k in eachindex(leftOrdinaryBoxes)
            criticalBoxNumber = incidences[i].leftAssignments[k]

            if criticalBoxNumber != 0
                criticalBox = incidences[i].boxes[criticalBoxNumber]
                push!(
                    edges,
                    topologyOrderOneEdge(leftOrdinaryBoxes[k], criticalBox),
                )
            end
        end

        for k in eachindex(rightOrdinaryBoxes)
            criticalBoxNumber = incidences[i].rightAssignments[k]

            if criticalBoxNumber != 0
                criticalBox = incidences[i].boxes[criticalBoxNumber]
                push!(
                    edges,
                    topologyOrderOneEdge(criticalBox, rightOrdinaryBoxes[k]),
                )
            end
        end
    end

    topologyAddVerticalEdges!(edges, V, criticalPoints, boxGroups)

    sort!(edges, lt = topologyEdgeIsLess)
    return edges
end


# A small example for interactive use.  This function is not run when the file
# is included.
function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    P = (x^2 + y^2 - 1) * x

    edges = connectCriticalBoxes(
        P;
        boxPrecision = 32,
        intervalPrecision = 16,
    )

    println(edges)
    println(length(edges))
end
