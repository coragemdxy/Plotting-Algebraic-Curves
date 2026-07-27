#=
Real isolation for the special points of a plane algebraic curve.

For a single polynomial P(x, y), the set P = 0 is normally infinite, so it
cannot be returned as finitely many boxes containing one root each.  This file
returns all the points on the special vertical fibers:

  1. all the points of H = 0 over a critical x-coordinate;
  2. all the intersections of a vertical component x = alpha with H = 0.

Each output box is ((x_left, x_right), (y_left, y_right)), and all boxes are
returned in one Set.
=#

using Nemo

# Make sure the polynomial is a nonzero QQ polynomial.
function requireBivariateQQPolynomial(P)
    R = parent(P)
    ngens(R) == 2 || throw(ArgumentError("P must have exactly two variables"))
    base_ring(R) == QQ || throw(ArgumentError("P must belong to QQ[x,y]"))
    iszero(P) && throw(ArgumentError("P must not be the zero polynomial"))
    return R
end

# Squarefree part of a multivariate polynomial over a characteristic-zero field.
function bivariateSquarefreePart(P)
    requireBivariateQQPolynomial(P)

    R = parent(P)
    n = ngens(R)
    G = P
    for i in 1:n
        dP = derivative(P, i)
        G = gcd(G, dP)
    end

    return divexact(P, G)
end

# Convert a QQ[x,y] polynomial depending only on variable `keep` to QQ[t].
function toUnivariate(P, keep::Int)
    keep in (1, 2) || throw(ArgumentError("keep must be 1 or 2"))
    S, t = polynomial_ring(QQ, keep == 1 ? "x" : "y")
    q = zero(S)

    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        other = keep == 1 ? e[2] : e[1]
        iszero(other) || return nothing
        q += c * t^e[keep]
    end
    return q
end

# Embed a QQ[t] polynomial in QQ[x,y] as a polynomial in variable `keep`.
function fromUnivariate(q, R, keep::Int)
    xs = gens(R)
    p = zero(R)
    for i in 0:degree(q)
        p += coeff(q, i) * xs[keep]^i
    end
    return p
end

# Evaluate a QQ[x,y] polynomial exactly at two algebraic numbers.
function evaluateAtAlgebraicPoint(P, a, b)
    Qb = parent(a)
    value = zero(Qb)
    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        value += Qb(c) * a^e[1] * b^e[2]
    end
    return value
end


# Vertical decomposition: P(x,y) = V(x) H(x,y)

# Return the content V(x) of P with respect to y, and the primitive part H.
function verticalDecomposition(P)
    R = requireBivariateQQPolynomial(P)
    Rx, t = polynomial_ring(QQ, "x")
    coefficients_in_y = Dict{Int, elem_type(Rx)}()

    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        coefficients_in_y[e[2]] =
            get(coefficients_in_y, e[2], zero(Rx)) + c * t^e[1]
    end

    nonzero_coefficients = filter(!iszero, collect(values(coefficients_in_y)))
    isempty(nonzero_coefficients) && error("internal error: no coefficient found")

    V = first(nonzero_coefficients)
    n = length(nonzero_coefficients)
    for a in 2:n
        V = gcd(V, nonzero_coefficients[a])
    end

    # A monic V gives a canonical decomposition and avoids irrelevant constants.
    V = divexact(V, leading_coefficient(V))
    Vxy = fromUnivariate(V, R, 1)
    H = divexact(P, Vxy)
    return V, H
end

# Get the leading coefficient of P for variable y as a polynomial in QQ[x].
function leadingCoefficientInY(P)
    S, t = polynomial_ring(QQ, "x")
    d = degree(P, 2)
    q = zero(S)

    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        if e[2] == d
            q += c * t^e[1]
        end
    end
    return q
end

# Make a univariate polynomial squarefree.
function univariateSquarefreePart(P)
    degree(P) <= 0 && return P
    return divexact(P, gcd(P, derivative(P)))
end


# Buchberger algorithm for QQ[x,y]

# True when exponent vector a is larger than b in the requested lex order.
function lexGreater(a, b, priority)
    for i in priority
        a[i] == b[i] || return a[i] > b[i]
    end
    return false
end

# Return coefficient and exponent vector of the lexicographic leading term.
function lexLeadingTerm(P, priority)
    iszero(P) && throw(ArgumentError("zero polynomial has no leading term"))
    cs = collect(coefficients(P))
    es = collect(exponent_vectors(P))
    best = 1
    for i in 2:length(es)
        if lexGreater(es[i], es[best], priority)
            best = i
        end
    end
    return cs[best], es[best]
end

function monomialFromExponent(R, e)
    xs = gens(R)
    return xs[1]^e[1] * xs[2]^e[2]
end

function monomialDivides(a, b)
    for i in eachindex(a)
        if a[i] > b[i]
            return false
        end
    end
    return true
end

# Multivariate division remainder using the requested lexicographic order.
function normalForm(P, basis, priority)
    R = parent(P)
    p = P
    r = zero(R)

    while !iszero(p)
        cp, ep = lexLeadingTerm(p, priority)
        reduced = false

        for g in basis
            iszero(g) && continue
            cg, eg = lexLeadingTerm(g, priority)
            if monomialDivides(eg, ep)
                multiplier = (cp / cg) * monomialFromExponent(R, ep .- eg)
                p -= multiplier * g
                reduced = true
                break
            end
        end

        if !reduced
            leading_term = cp * monomialFromExponent(R, ep)
            r += leading_term
            p -= leading_term
        end
    end
    return r
end

function sPolynomial(F, G, priority)
    R = parent(F)
    cf, ef = lexLeadingTerm(F, priority)
    cg, eg = lexLeadingTerm(G, priority)
    lcm_exp = max.(ef, eg)
    return monomialFromExponent(R, lcm_exp .- ef) * F / cf -
           monomialFromExponent(R, lcm_exp .- eg) * G / cg
end

function makeMonic(P, priority)
    iszero(P) && return P
    c, _ = lexLeadingTerm(P, priority)
    return P / c
end

# Compute a lexicographic Groebner basis by Buchberger's algorithm.
function buchbergerBasis(polynomials; priority = [1, 2])
    # Algorithm 15 starts with G := F, so keep the input polynomials unchanged.
    basis = collect(polynomials)
    isempty(basis) && return basis

    pairs = [(i, j) for i in 1:length(basis) for j in (i+1):length(basis)]
    cursor = 1
    while cursor <= length(pairs)
        i, j = pairs[cursor]
        cursor += 1
        r = normalForm(sPolynomial(basis[i], basis[j], priority), basis, priority)
        if !iszero(r)
            old_length = length(basis)
            push!(basis, r)
            for k in 1:old_length
                push!(pairs, (k, old_length + 1))
            end
        end
    end

    # Return G reordered by the requested monomial order (largest LM first).
    sort!(basis; lt = (f, g) -> begin
        _, ef = lexLeadingTerm(f, priority)
        _, eg = lexLeadingTerm(g, priority)
        lexGreater(ef, eg, priority)
    end)
    return basis
end

# Get an eliminant in x or y from a lex Groebner basis.
function eliminationPolynomial(F, G, keep::Int)
    # To keep x, eliminate y first (priority y > x), and conversely for y.
    priority = keep == 1 ? [2, 1] : [1, 2]
    basis = buchbergerBasis([F, G]; priority = priority)
    candidates = []
    for g in basis
        q = toUnivariate(g, keep)
        if q !== nothing && q != 0
            push!(candidates, q)
        end
    end

    if isempty(candidates)
        # A zero-dimensional elimination ideal must contain a nonzero
        # univariate polynomial.
        throw(ArgumentError("the system is not zero-dimensional"))
    end

    sort!(candidates, by = degree)
    q = first(candidates)
    return divexact(q, gcd(q, derivative(q)))
end


# Exact algebraic solutions and rational isolating boxes

# Isolate one real QQBar number in an interval of width at most 2^-precise.
function isolateAlgebraicNumber(a, precise::Int)
    precise >= 0 || throw(ArgumentError("precise must be nonnegative"))
    Qb = parent(a)
    bound = ZZ(1)
    while !(Qb(-bound) <= a <= Qb(bound))
        bound *= 2
    end

    left = QQ(-bound)
    right = QQ(bound)
    target = QQ(1, ZZ(2)^precise)

    while right - left > target
        middle = (left + right) / 2
        if a < Qb(middle)
            right = middle
        elseif Qb(middle) < a
            left = middle
        else
            return (middle, middle)
        end
    end
    return (left, right)
end

# Return all distinct exact real solutions of F = G = 0.
function zeroDimensionalPoints(F, G, precise::Int)
    parent(F) == parent(G) || throw(ArgumentError("F and G must have the same parent"))
    iszero(F) && throw(ArgumentError("F must not be zero"))
    iszero(G) && throw(ArgumentError("G must not be zero"))

    tx = eliminationPolynomial(F, G, 1)
    ty = eliminationPolynomial(F, G, 2)
    Qb = algebraic_closure(QQ)
    xs = unique(filter(is_real, roots(Qb, tx)))
    ys = unique(filter(is_real, roots(Qb, ty)))
    points = NamedTuple[]

    # The eliminants give coordinate candidates.  Exact substitution pairs the
    # correct x and y values and rejects the Cartesian-product false positives.
    for a in xs, b in ys
        if iszero(evaluateAtAlgebraicPoint(F, a, b)) &&
           iszero(evaluateAtAlgebraicPoint(G, a, b))
            box = (
                isolateAlgebraicNumber(a, precise),
                isolateAlgebraicNumber(b, precise),
            )
            push!(points, (x = a, y = b, box = box))
        end
    end
    sort!(points, by = p -> (Float64(p.x), Float64(p.y)))
    return points
end

samePoint(p, q) = p.x == q.x && p.y == q.y


# Curve-isolation algorithm

function realRootsAsAlgebraic(q)
    degree(q) <= 0 && return QQBarFieldElem[]
    Qb = algebraic_closure(QQ)
    squarefree_q = divexact(q, gcd(q, derivative(q)))
    return unique(filter(is_real, roots(Qb, squarefree_q)))
end

# Isolate the finite topology data associated with P(x,y) = 0.
function bivariateRealIsolation(P, precise::Int=32)
    requireBivariateQQPolynomial(P)
    precise >= 0 || throw(ArgumentError("precise must be nonnegative"))
    Psf = bivariateSquarefreePart(P)
    V, H = verticalDecomposition(Psf)

    # C contains the x-coordinates of the critical points of H = 0.
    # The system H = C = 0 returns every point on these critical fibers,
    # including the points which do not satisfy derivative(H, 2) = 0.
    special_fiber_points = NamedTuple[]
    if degree(H, 2) > 0
        Hy = derivative(H, 2)
        if !iszero(Hy)
            C = eliminationPolynomial(H, Hy, 1)
            L = leadingCoefficientInY(H)
            E = univariateSquarefreePart(C * L)

            if degree(E) > 0
                R = parent(Psf)
                Exy = fromUnivariate(E, R, 1)
                special_fiber_points = zeroDimensionalPoints(H, Exy, precise)
            end
        end
    end

    vertical_components = NamedTuple[]
    all_intersections = NamedTuple[]
    vertical_roots = realRootsAsAlgebraic(V)

    if !isempty(vertical_roots) && degree(H, 2) > 0
        R = parent(Psf)
        Vxy = fromUnivariate(V, R, 1)
        all_intersections = zeroDimensionalPoints(Vxy, H, precise)
    end

    for alpha in vertical_roots
        x_interval = isolateAlgebraicNumber(alpha, precise)
        points_on_component = [p for p in all_intersections if p.x == alpha]
        sort!(points_on_component, by = p -> Float64(p.y))
        push!(vertical_components, (
            x_interval = x_interval,
            intersections = [p.box for p in points_on_component],
        ))
    end

    # A point can belong to a special fiber and a vertical component.
    # Keep it only once in the final result.
    special_fiber_points = [
        p for p in special_fiber_points
        if !any(q -> samePoint(p, q), all_intersections)
    ]

    boxes = [p.box for p in special_fiber_points]
    for component in vertical_components
        append!(boxes, component.intersections)
    end
    sort!(boxes, by = box -> (
        (box[1][1] + box[1][2]) / 2,
        (box[2][1] + box[2][2]) / 2,
    ))
    return Set(boxes)
end


function main()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 1) * (x^2 + y^2 - 4) * x
    boxes = bivariateRealIsolation(f, 16)
    println("f = ", f)
    display(boxes)
end
