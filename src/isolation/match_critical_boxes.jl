#=
Match the boxes returned by `bivariateRealIsolation` with the critical
x-coordinates returned by `getPointsCriticalWithIntervals`.

The output is a sorted vector.  Each element has the form

    ((alpha, (x_left, x_right)), (y_left, y_right))

where `alpha` is the exact critical x-coordinate, `(x_left, x_right)` is its
rational isolating interval, and `(y_left, y_right)` is the y-interval from the
original bivariate box.
=#

using Nemo

# Load the two producers when this file is included on its own.  The guards
# avoid redefining their functions when the caller has already loaded them.
if !isdefined(@__MODULE__, :bivariateRealIsolation)
    include(joinpath(@__DIR__, "bivariate_real_isolation.jl"))
end
if !isdefined(@__MODULE__, :getPointsCritical) ||
   !isdefined(@__MODULE__, :getPointsCriticalWithIntervals)
    include(joinpath(@__DIR__, "qbar_critical_points.jl"))
end


# Check that an interval has its endpoints in increasing order.
function requireOrderedInterval(interval, name::AbstractString)
    length(interval) == 2 ||
        throw(ArgumentError("$name must contain exactly two endpoints"))

    left, right = interval
    left <= right ||
        throw(ArgumentError("$name has endpoints in the wrong order: $interval"))

    return interval
end


# True when `inner` is contained in `outer`, including their boundaries.
function intervalContains(outer, inner)
    requireOrderedInterval(outer, "outer interval")
    requireOrderedInterval(inner, "inner interval")

    outerLeft, outerRight = outer
    innerLeft, innerRight = inner
    return outerLeft <= innerLeft && innerRight <= outerRight
end





#Match box and interval and replace the x-interval of every bivariate box by the corresponding exact critical point and its isolating interval.

function matchBoxesToCriticalIntervals(boxes, criticalIntervals)
    criticalData = collect(criticalIntervals)
    sort!(criticalData, by = item -> item[1])

    for (i, item) in enumerate(criticalData)
        length(item) == 2 ||
            throw(ArgumentError("critical entry $i must have the form (alpha, interval)"))
        requireOrderedInterval(item[2], "critical interval $i")
    end

    # The intervals produced by getPointsCriticalWithIntervals are disjoint.
    # Check this explicitly so that every bivariate box has at most one match.
    for i in 1:(length(criticalData) - 1)
        currentRight = criticalData[i][2][2]
        nextLeft = criticalData[i + 1][2][1]
        currentRight < nextLeft ||
            throw(ArgumentError(
                "critical intervals $i and $(i + 1) overlap or touch",
            ))
    end

    result = []

    for (i, box) in enumerate(boxes)
        length(box) == 2 ||
            throw(ArgumentError("box $i must have the form (x_interval, y_interval)"))

        xInterval = requireOrderedInterval(box[1], "x-interval of box $i")
        yInterval = requireOrderedInterval(box[2], "y-interval of box $i")
        matchingIndices = Int[]

        for (j, (_, criticalInterval)) in enumerate(criticalData)
            if intervalContains(criticalInterval, xInterval)
                push!(matchingIndices, j)
            end
        end

        if isempty(matchingIndices)
            throw(ArgumentError(
                "the x-interval $xInterval of box $i is not contained in any critical interval",
            ))
        elseif length(matchingIndices) > 1
            throw(ArgumentError(
                "the x-interval $xInterval of box $i matches more than one critical interval",
            ))
        end

        criticalPoint, criticalInterval = criticalData[only(matchingIndices)]
        push!(result, ((criticalPoint, criticalInterval), yInterval))
    end

    sort!(result, lt = (firstBox, secondBox) -> begin
        firstPoint = firstBox[1][1]
        secondPoint = secondBox[1][1]

        if firstPoint != secondPoint
            return firstPoint < secondPoint
        end

        firstY = firstBox[2]
        secondY = secondBox[2]
        if firstY[1] != secondY[1]
            return firstY[1] < secondY[1]
        end
        return firstY[2] < secondY[2]
    end)

    return result
end



#Compute the bivariate boxes and critical x-coordinates of `P`, match them, and return boxes in the form `((alpha, critical_x_interval), original_y_interval)`.
function criticalPointBoxes(P; boxPrecision::Int = 64,intervalPrecision::Int = 32)
    boxPrecision >= 0 ||
        throw(ArgumentError("boxPrecision must be nonnegative"))
    intervalPrecision >= 0 ||
        throw(ArgumentError("intervalPrecision must be nonnegative"))

    isdefined(@__MODULE__, :bivariateRealIsolation) ||
        throw(ArgumentError("load bivariate_real_isolation.jl first"))
    isdefined(@__MODULE__, :getPointsCritical) ||
        throw(ArgumentError("load Qbar_critical_points.jl first"))
    isdefined(@__MODULE__, :getPointsCriticalWithIntervals) ||
        throw(ArgumentError("load Qbar_critical_points.jl first"))

    boxes = bivariateRealIsolation(P, boxPrecision)
    criticalPoints = getPointsCritical(P)
    criticalIntervals =
        getPointsCriticalWithIntervals(criticalPoints, intervalPrecision)

    return matchBoxesToCriticalIntervals(boxes, criticalIntervals)
end
