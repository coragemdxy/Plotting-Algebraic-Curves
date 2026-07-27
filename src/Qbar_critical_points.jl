#=
This program is used to get all the roots in Qbar of a polynomial QQ.
=#

using Nemo

function coeffInVar(F,a,b)
    R = parent(F)
    xs = gens(R)
    n = length(xs)
    d = degree(F, 2)
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

function convertToQQ(P)
    if iszero(P)
        return zero(QQ)
    end
    return collect(coefficients(P))[1]
end

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

function sylvesterResultant(F, G, a)
    S = sylvesterMatrix(F, G, a)
    Q =  det(S)
    T,x = polynomial_ring(QQ,"x") 
    f = zero(T)
    for (g,e) in zip(coefficients(Q), exponent_vectors(Q))
            f += g*x^(e[1])
    end
    return f
end

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

function getLeadingCoef(P)
    R = parent(P)
    d = degree(P, 2)
    
    S,x = polynomial_ring(QQ,"x") 
    f = zero(S)
    for (g,e) in zip(coefficients(P), exponent_vectors(P))
        if e[2] == d
            f += g*x^(e[1])
        end
    end
    return f
end

function getRootsOfLeadingCoef(P)
    Qb = algebraic_closure(QQ)
    Q = getLeadingCoef(P)
    s = roots(Qb,Q)
    s1 = []
    for a in s
        if is_real(a)
            push!(s1,a)
        end
    end
    return s1
end

function getRootsOfDiscriminant(P)
    Qb = algebraic_closure(QQ)
    Q = derivative(P,2)
    Q1 = sylvesterResultant(P,Q,2)
    s = roots(Qb,Q1)
    s2 = []
    for a in s
        if is_real(a)
            push!(s2,a)
        end
    end
    return s2
end

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
    set3 = unique(union(set1, set2))
    sort!(set3)
    return set3
end

function getPointsCriticalWithIntervals(points, precise)
    isempty(points) && return []

    precise >= 0 ||
        throw(ArgumentError("precise must be nonnegative"))

    points1 = sort(unique(points))
    Qb = parent(points1[1])
    result = []

    for a in points1
        bound = ZZ(1)

        while !(Qb(-bound) < a < Qb(bound))
            bound *= 2
        end

        left = QQ(-bound)
        right = QQ(bound)
        push!(result, (a, (left, right)))
    end

    targetWidth = QQ(1, ZZ(2)^precise)

    while true
        separated = true
        sufficientlyPrecise = true

        if length(result) >= 2
            for i in 1:(length(result) - 1)
                right1 = result[i][2][2]
                left2 = result[i + 1][2][1]

                if right1 >= left2
                    separated = false
                    break
                end
            end
        end

        for (_, interval) in result
            left, right = interval

            if right - left > targetWidth
                sufficientlyPrecise = false
                break
            end
        end

        if separated && sufficientlyPrecise
            return result
        end

        result1 = []

        for (a, interval) in result
            left, right = interval
            middle = (left + right) / 2

            if a < Qb(middle)
                right = middle
            elseif Qb(middle) < a
                left = middle
            else
                # 根恰好等于有理中点时，构造一个仍严格包含它的小区间
                width = (right - left) / 4
                left = middle - width
                right = middle + width
            end

            push!(result1, (a, (left, right)))
        end

        result = result1
    end
end

function main()
    R,(x,y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 3)
    println(f)
    s = getPointsCritical(f)
    s1 = getPointsCriticalWithIntervals(s,32)
    println(s1)
end
