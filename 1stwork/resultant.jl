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