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
    println(f)
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
    return det(S)
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
    xs = gens(R)
    x = xs[1]
    d = degree(P, 2)
    f = 0
    for (g,e) in zip(coefficients(P), exponent_vectors(P))
        if e[2] == d
            f += g*x^(e[1])
        end
    end
    println(f)
    return f
end

function getRootsOfLeadingCoef(P)
    R = parent(P)
    Q = getLeadingCoef(P)
    fac1 = factor(Q)
    println(fac1)
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

function getRootsOfDiscriminant(P)
    R = parent(P)
    Q = derivative(P,2)
    Q1 = sylvesterResultant(P,Q,2)
    Q2 = factor(Q1)
    println(Q2)
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

function getSamplePoints(P)
    if iszero(P)
        throw(ArgumentError("P can not be a zero polynomial."))
    end

    R = parent(P)        
    if base_ring(R) != QQ
        throw(ArgumentError("P is not a rationnal polynomial"))
    end

    P1 = squareFree(P)
    println(P1)
    set1 = getRootsOfLeadingCoef(P1)

    set2 = getRootsOfDiscriminant(P1)

    return union(set1, set2)
end
    
function main1()
    R,(x,y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*x
    println(f)
    s = getSamplePoints(f)
    println(s)
end