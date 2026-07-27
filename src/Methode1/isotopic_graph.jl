#=
This program is used to connect sample points of a polynomial, and make sure that the graph connected is 
isotopic to the graph related to the polynomial.
=#

using Nemo

##Get the coefficients of variation a, degree b
function coeffInVar(F,a,b)
    R = parent(F)
    xs = gens(R)
    n = length(xs)
    f = zero(R)
    for (g,e) in zip(coefficients(F), exponent_vectors(F))
        if e[a] == b
            term = R(g)
            for i in 1:n
                if i != a
                    term *= xs[i]^e[i]
                end
            end
            f += term
        end
    end
    return f
end

#convert a float to QQ
function convertToQQ(P)
    if iszero(P)
        return zero(QQ)
    end
    return collect(coefficients(P))[1]
end

#Calculer the sylvesterMatrix of F,G for the variable a
function sylvesterMatrix(F, G, a)
    parent(F) == parent(G) ||
        throw(ArgumentError("$F and $G must belong to the same polynomial ring"))

    iszero(F) && throw(ArgumentError("$F is a zero polynomial"))
    iszero(G) && throw(ArgumentError("$G is a zero polynomial"))

    m = degree(F,a)
    n = degree(G,a)
    R = parent(F)
    N = m + n

    S = zero_matrix(R, N, N)
    for row in 1:n
        for j in 0:m
            S[row, row + j] = coeffInVar(F,a,m - j)
        end
    end
    for row in 1:m
        for j in 0:n
            S[n + row, row + j] = coeffInVar(G,a,n - j)
        end
    end

    return S
end

#Calculer the resultant of F,G for the variable a
function sylvesterResultant(F, G, a)
    S = sylvesterMatrix(F, G, a)
    return det(S)
end

#Make a polynomial P squarefree
function squareFree(P)
    if iszero(P)
        throw(ArgumentError("zero polynomial has no squarefree part"))
    end

    R = parent(P)
    n = ngens(R)
    G = P
    for i in 1:n
        dP = derivative(P, i)
        G = gcd(G, dP)
    end

    return divexact(P, G)
end

function squareFreeUnivar(P)
    if iszero(P)
        throw(ArgumentError("zero polynomial has no squarefree part"))
    end

    dP = derivative(P)
    G = gcd(P,dP)
    return divexact(P, G)
end

#Get the the LeadingCoef for variable y
function getLeadingCoef(P)
    R = parent(P)
    xs = gens(R)
    x = xs[1]
    d = degree(P, 2)
    f = 0
    for (g,e) in zip(coefficients(P), exponent_vectors(P))
        if e[2] == d
            f += g*x^(e[1])
        end
    end
    return f
end

#Get the roots of Leading coefficients of variable y
function getRootsOfLeadingCoef(P)
    R = parent(P)
    Q = getLeadingCoef(P)
    fac1 = factor(Q)
    s1 = Set{elem_type(QQ)}()
    for (g,e) in fac1
        if degree(g,1) >=2
            throw(error("There exists an irrational roots"))
        elseif degree(g,1) == 1
            a = -coeffInVar(g,1,0)
            a0 = convertToQQ(a)
            a= coeffInVar(g,1,1)
            a1 = convertToQQ(a)
            r = a0/a1
            push!(s1, r)
        end
    end
    return s1 
end

#Get the roots of Discriminant of P for variable y
function getRootsOfDiscriminant(P)
    R = parent(P)
    Q = derivative(P,2)
    Q1 = sylvesterResultant(P,Q,2)
    Q2 = factor(Q1)
    s2 = Set{elem_type(QQ)}()
    for (g,e) in Q2
        if degree(g,1) >=2
            throw(error("There exists an irrational roots"))
        elseif degree(g,1) == 1
            a = -coeffInVar(g,1,0)
            a0 = convertToQQ(a)
            a= coeffInVar(g,1,1)
            a1 = convertToQQ(a)
            r = a0/a1
            push!(s2, r)
        end
    end
    return s2

end

#Get all the x of sample points.
function getPointsCritical(P)
    if iszero(P)
        throw(ArgumentError("P can not be a zero polynomial."))
    end

    R = parent(P)        
    if base_ring(R) != QQ
        throw(ArgumentError("P is not a rationnal polynomial"))
    end

    P1 = squareFree(P)
    set1 = getRootsOfLeadingCoef(P1)

    set2 = getRootsOfDiscriminant(P1)

    return union(set1, set2)
end

#Calculer the bound of roots of a univariate polynomial
function boundOfRoots(P)
    array = calSequenceOfCoef(P)
    n = length(array)
    S = zero(QQ)
    for i in 1:(n-1)
        S += abs(array[i]/array[n])
    end
    M = max(one(QQ), S)

    mx = zero(QQ)
    for i in 1:(n-1)
        if mx < abs(array[i]/array[n])
            mx = abs(array[i]/array[n])
        end
    end
    N = one(QQ) + mx

    return min(M, N)
end

#Calculer the number of sign variation of a sequence
function signVariation(list)
    array = filter(x -> x != 0,list)
    n = length(array)

    if n <= 1
        return 0
    end
    num = 0
    for i in 1:(n-1)
        if array[i]*array[i+1] < 0
            num += 1
        end
    end
    return num
end

#Calculer the sequence of coefficients of a polynomial
function calSequenceOfCoef(P)
    array = []
    n = degree(P)

    for i in 0:n
        push!(array,coeff(P,i))
    end
    return array
end

#Make mobius transformation for a polynomial and an interval
function mobiusTransformation(P, interval)
    R = parent(P)
    x = gens(R)[1]

    arraycoef = calSequenceOfCoef(P)
    n = length(arraycoef)
    a = interval[1]
    b = interval[2]

    Q = zero(QQ)
    for i in 0:(n-1)
        Q += arraycoef[i+1]*(a+b*x)^(i)*(1+x)^(n-i-1)
    end
    return Q
end

#Use Descartes’ law of signs adn Mobius transformation to calculer the interval of every root of a polynomial
function realIsolationPart(P, interval,precise,times)
    res = []

    precise1 =  one(QQ) / (QQ(2)^precise)
    maxtimes = 500
    if times > maxtimes
        throw(error("The times of recurrence is over the maximum"))
    end
    Q = mobiusTransformation(P,interval)
    num = signVariation(calSequenceOfCoef(Q))
    if num == 1 && (interval[2]-interval[1]<precise1)
        push!(res,interval)
    elseif num == 0
        return res
    else
        mid = (interval[1]+interval[2])/2
        if P(mid) == 0
            push!(res, [mid,mid])
        end
        interval1 = [interval[1],mid]
        interval2 = [mid,interval[2]]
        append!(res, realIsolationPart(P, interval1,precise, times+1))
        append!(res, realIsolationPart(P, interval2,precise,times+1))

    end
    return res
end

#Get the root of component x-a
function realIsolationVertical(P,a,precise)
    R,(x,y) = polynomial_ring(QQ,["x","y"])
    H = divexact(P, x - QQ(a))
    Ry, x = polynomial_ring(QQ, "y")
    Q = evaluate(H, [QQ(a), x])
    res = []
    if degree(Q) <=0
        res = realIsolationVertical(P,a,precise)
    else
        p = squareFreeUnivar(Q)
        M = boundOfRoots(p)

        interval = [-M-1, M+1]
        push!(res,(a,[-M-1,-M-1]))
        res = realIsolationPart(p,interval,precise,0)
        push!(res,(a,[M+1,M+1]))
    end
    return res
end


#Get the interval of roots of P at x=a, by using real isolation with the precise.
function getRealIsolationOfRoots(P,a,precise)
    Ry, x = polynomial_ring(QQ, "y")
    Q = evaluate(P, [QQ(a), x])

    res = []
    if degree(Q) <=0
        res = realIsolationVertical(P,a,precise)
    else
        p = squareFreeUnivar(Q)
        M = boundOfRoots(p)
        interval = [-M-1, M+1]
        res =realIsolationPart(p,interval,precise,0)
    end
    return res

end

#Get all the sample points of P with the precise
function getSamplePoints(P,precise)
    array = getPointsCritical(P)

    array1 = sort!(collect(array))
    n = length(array1)
    array2 = [array1[1]-1]
    for i in 1:(n-1)
        push!(array2,array1[i])
        push!(array2, (array1[i]+array1[i+1])/2)
    end
    push!(array2, array1[n])
    push!(array2,array1[n]+1)

    array3 = []
    for a in array2
        array4 = getRealIsolationOfRoots(P,a,precise)
        for i in array4
            push!(array3,(a,i))
        end
    end

    return array3
end

#Get the middle point of an interval
function getMiddleOfInterval(interval)
    return (interval[1] + interval[2]) / 2
end

#Get the y coordinate of a sample point
#A sample point has the form (x, interval)
function getYOfSamplePoint(point)
    return getMiddleOfInterval(point[2])
end

#Get all x coordinates from sample points
function getXCoordinatesOfSamplePoints(points)
    array1 = []
    for i in points
        if !(i[1] in array1)
            push!(array1, i[1])
        end
    end
    sort!(array1)
    return array1
end

#Get the indices of points whose x coordinate is a
function getIndicesOfPointsAtX(points, a)
    array1 = []
    for i in eachindex(points)
        if points[i][1] == a
            push!(array1, i)
        end
    end

    sort!(array1, by = i -> getYOfSamplePoint(points[i]))
    return array1
end

#Judge whether a is a critical x coordinate
function isCriticalPoint(a, criticals)
    return a in criticals
end

#Judge whether x = a is a vertical component of P = 0
function isVerticalComponentAtX(P, a)
    Ry, x = polynomial_ring(QQ, "y")
    Q = evaluate(P, [QQ(a), x])
    if degree(Q) <=0
        return true
    end
    return false
end

#Get sorted root intervals of P(a,y)=0
function getSortedRootsAtX(P, a, precise)
    array1 = getRealIsolationOfRoots(P, a, precise)
    sort!(array1, by = i -> getMiddleOfInterval(i))
    return array1
end

#Find the nearest critical point according to the y coordinate
function getNearestPoint(points, array1, y)
    if length(array1) == 0
        throw(error("There is no point in array1."))
    end

    index = array1[1]
    distance = abs(getYOfSamplePoint(points[index]) - y)

    for i in array1
        distance1 = abs(getYOfSamplePoint(points[i]) - y)
        if distance1 < distance
            distance = distance1
            index = i
        end
    end

    return index
end

#Add an edge to the adjacency matrix
function addEdge!(M, i, j)
    M[i,j] = 1
    M[j,i] = 1
end

#Connect one ordinary vertical line with one critical vertical line
function connectOrdinaryAndCritical!(M, P, points, ordinaryX, ordinaryIndices, criticalX, criticalIndices, precise)
    auxX = (ordinaryX + criticalX) / 2

    array1 = getSortedRootsAtX(P, auxX, precise)

    if length(array1) != length(ordinaryIndices)
        throw(error("The number of roots at auxX is different from the number of roots at ordinaryX."))
    end

    for i in eachindex(ordinaryIndices)
        ordinaryIndex = ordinaryIndices[i]

        y = getMiddleOfInterval(array1[i])

        criticalIndex = getNearestPoint(points, criticalIndices, y)

        addEdge!(M, ordinaryIndex, criticalIndex)
    end
end

#Connect points on the same vertical component
function connectVerticalComponent!(M, points, indices)
    if length(indices) <= 1
        return
    end

    sort!(indices, by = i -> getYOfSamplePoint(points[i]))

    for k in 1:(length(indices)-1)
        addEdge!(M, indices[k], indices[k+1])
    end
end

#Get the adjacency matrix of the graph
function getAdjacencyMatrix(P, points, precise)
    array1 = sort!(collect(getPointsCritical(P)))
    array2 = getXCoordinatesOfSamplePoints(points)

    n = length(points)
    M = zeros(Int, n, n)

    for i in 1:(length(array2)-1)
        x1 = array2[i]
        x2 = array2[i+1]

        array3 = getIndicesOfPointsAtX(points, x1)
        array4 = getIndicesOfPointsAtX(points, x2)

        b1 = isCriticalPoint(x1, array1)
        b2 = isCriticalPoint(x2, array1)

        if b1 && !b2
            connectOrdinaryAndCritical!(M, P, points, x2, array4, x1, array3, precise)
        elseif !b1 && b2
            connectOrdinaryAndCritical!(M, P, points, x1, array3, x2, array4, precise)
        end
    end

    for x in array2
        if isCriticalPoint(x, array1) && isVerticalComponentAtX(P, x)
            indices = getIndicesOfPointsAtX(points, x)
            connectVerticalComponent!(M, points, indices)
        end
    end

    return M
end

#Convert adjacency matrix to the list of edges
function getEdges(M)
    array1 = []
    n = size(M, 1)

    for i in 1:n
        for j in (i+1):n
            if M[i,j] == 1
                push!(array1, (i,j))
            end
        end
    end

    return array1
end

function main()
    R,(x,y) = polynomial_ring(QQ,["x","y"])

    precise = 60

    f = y^2 + x^2-1

    points, M, edges = getGraph(f, precise)

    println("points = ")
    println(points)

    println("adjacency matrix = ")
    println(M)

    println("edges = ")
    println(edges)
end