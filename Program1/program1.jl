using Nemo

function sylvester_matrix(F, G)
    parent(F) == parent(G) ||
        throw(ArgumentError("$F and $G must belong to the same polynomial ring"))

    iszero(F) && throw(ArgumentError("$F is a zero polynomial"))
    iszero(G) && throw(ArgumentError("$G is a zero polynomial"))

    m = degree(F)
    n = degree(G)
    R = base_ring(parent(F))
    N = m + n

    S = zero_matrix(R, N, N)
    for row in 1:n
        for j in 0:m
            S[row, row + j] = coeff(F, m - j)
        end
    end
    for row in 1:m
        for j in 0:n
            S[n + row, row + j] = coeff(G, n - j)
        end
    end

    return S
end

function sylvester_resultant(F, G)
    S = sylvester_matrix(F, G)
    return det(S)
end

function getRootsOfLeadingCoef(P)
    Q = leading_coefficient(P)
    fac1 = factor(Q)
    s1 = Set{elem_type(QQ)}()
    for a in fac1
        if degree(a) >=2
            throw(error("There exists a roots irrational"))
        end
        push!(s1, roots(a))
    end
    return s1 
end

function getRootsOfDiscriminant(P)
    Q = derivative(P)
    Q1 = sylvester_resultant(P,Q)
    Q2 = factor(Q1)
    println(Q2)
    s2 = Set{elem_type(QQ)}()
    for (g,e) in Q2
        if degree(g) >=2
            throw(error("There exists a roots irrational"))
        elseif degree(g) == 1
            r = -coeff(g, 0) / coeff(g, 1)
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
    Q = base_ring(R)

    if !(Q isa PolyRing)
        throw(ArgumentError("P must be in R[y], where R is a polynomial ring."))
    end
    if base_ring(Q) != QQ
        throw(ArgumentError("The coefficient ring must be QQ[x]."))
    end

    set1 = getRootsOfLeadingCoef(P)

    set2 = getRootsOfDiscriminant(P)

    return union(set1, set2)
end
    
function main1()
    Q,x = polynomial_ring(QQ, "x")
    R,y  = polynomial_ring(Q, "y")
    f = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)
    s = getSamplePoints(f)
    println(s)
end