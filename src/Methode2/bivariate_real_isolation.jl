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
import AlgebraicSolving

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

# Rewrite a polynomial so that the variable to eliminate is the first variable
# and the variable to keep is the second one.  AlgebraicSolving eliminates the
# first variable block, while callers of eliminationPolynomial specify the
# variable which must be retained.
function polynomialForBlockElimination(P, S, eliminate::Int, keep::Int)
    eliminatedVariable, keptVariable = gens(S)
    q = zero(S)

    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        q += c *
            eliminatedVariable^e[eliminate] *
            keptVariable^e[keep]
    end
    return q
end

# Get an eliminant in x or y with AlgebraicSolving's F4 block-elimination
# algorithm.  Local resultants are used later only for matching coordinates
# inside one irreducible x-block.
function eliminationPolynomial(F, G, keep::Int)
    keep in (1, 2) || throw(ArgumentError("keep must be 1 or 2"))
    parent(F) == parent(G) ||
        throw(ArgumentError("F and G must have the same parent"))
    iszero(F) && throw(ArgumentError("F must not be zero"))
    iszero(G) && throw(ArgumentError("G must not be zero"))

    eliminate = keep == 1 ? 2 : 1
    S, _ = polynomial_ring(
        QQ,
        ["eliminated", "kept"];
        internal_ordering = :degrevlex,
    )
    reorderedF = polynomialForBlockElimination(F, S, eliminate, keep)
    reorderedG = polynomialForBlockElimination(G, S, eliminate, keep)
    ideal = AlgebraicSolving.Ideal([reorderedF, reorderedG])
    basis = AlgebraicSolving.eliminate(ideal, 1)

    T, t = polynomial_ring(QQ, keep == 1 ? "x" : "y")
    candidates = elem_type(T)[]
    for g in basis
        q = zero(T)
        for (c, e) in zip(coefficients(g), exponent_vectors(g))
            iszero(e[1]) ||
                error(
                    "internal error: block elimination retained " *
                    "the eliminated variable",
                )
            q += c * t^e[2]
        end
        iszero(q) || push!(candidates, q)
    end

    if isempty(candidates)
        throw(ArgumentError("the system is not zero-dimensional"))
    end

    sort!(candidates, by = degree)
    return univariateSquarefreePart(first(candidates))
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

# Return the distinct real roots of a univariate QQ polynomial.
function realRootsAsAlgebraic(q)
    degree(q) <= 0 && return QQBarFieldElem[]
    Qb = algebraic_closure(QQ)
    squarefree_q = univariateSquarefreePart(q)
    return unique(filter(is_real, roots(Qb, squarefree_q)))
end

# Evaluate a bivariate QQ polynomial by ball arithmetic.  A ball which excludes
# zero certifies that the corresponding Cartesian-product candidate is false.
function evaluateAtArbPoint(P, a, b)
    B = parent(a)
    value = zero(B)
    for (c, e) in zip(coefficients(P), exponent_vectors(P))
        value += B(c) * a^e[1] * b^e[2]
    end
    return value
end

function pairwiseDisjoint(balls)
    for i in 1:length(balls), j in (i + 1):length(balls)
        overlaps(balls[i], balls[j]) && return false
    end
    return true
end

# Use inexpensive certified interval rejection before constructing a separating
# projection.  True solutions are never removed: only balls on which H excludes
# zero are discarded.
function arbCandidatePairs(H, xs, ys)
    candidates = [(i, j) for i in eachindex(xs) for j in eachindex(ys)]
    for bits in (64, 128, 256, 512, 1024)
        B = ArbField(bits)
        xballs = B.(xs)
        yballs = B.(ys)
        candidates = [
            (i, j) for (i, j) in candidates
            if contains_zero(evaluateAtArbPoint(H, xballs[i], yballs[j]))
        ]
        isempty(candidates) && break
    end
    return candidates
end

# Rewrite lambda^d H(x, (t-x)/lambda) as a polynomial in x over QQ[t].
# Its resultant with q(x) is the exact projection polynomial for
# t = x + lambda*y.
function separatingProjectionData(H, q, lambda::Int)
    iszero(lambda) && throw(ArgumentError("lambda must not be zero"))

    T, t = polynomial_ring(QQ, "t")
    S, x = polynomial_ring(T, "x")
    qNested = zero(S)
    for i in 0:degree(q)
        qNested += T(coeff(q, i)) * x^i
    end

    d = degree(H, 2)
    transformed = zero(S)
    for (c, e) in zip(coefficients(H), exponent_vectors(H))
        transformed += T(c) *
            x^e[1] *
            (t - x)^e[2] *
            lambda^(d - e[2])
    end

    labelPolynomial = resultant(qNested, transformed)
    iszero(labelPolynomial) && return nothing
    labelPolynomial = univariateSquarefreePart(labelPolynomial)

    # A primitive pseudo-remainder sequence is much smaller than a Groebner
    # basis here.  When it reaches degree one it gives A(t)x + C(t) = 0.
    # Record every nonconstant content removed from a remainder: the relation
    # is used at a root tau only after certifying that all these contents are
    # nonzero at tau.
    firstPolynomial = transformed
    secondPolynomial = qNested
    if degree(firstPolynomial) < degree(secondPolynomial)
        firstPolynomial, secondPolynomial =
            secondPolynomial, firstPolynomial
    end

    guards = elem_type(T)[]
    while degree(secondPolynomial) > 1
        remainder = -pseudorem(firstPolynomial, secondPolynomial)
        iszero(remainder) && return nothing
        commonContent = content(remainder)
        if degree(commonContent) > 0
            push!(guards, commonContent)
            remainder = divexact(remainder, S(commonContent))
        end
        firstPolynomial, secondPolynomial =
            secondPolynomial, remainder
    end

    degree(secondPolynomial) == 1 || return nothing
    return (
        label = labelPolynomial,
        relation = secondPolynomial,
        guards = guards,
    )
end

# Certify a bijection between the retained (x,y) candidates and the real roots
# of a separating linear projection t=x+lambda*y.
function certifyCandidatePairs(H, q, xs, ys, candidates, lambda::Int)
    projection = separatingProjectionData(H, q, lambda)
    projection === nothing && return nothing

    ts = realRootsAsAlgebraic(projection.label)
    length(ts) == length(candidates) || return nothing

    relation = projection.relation
    coefficientOfX = coeff(relation, 1)
    constantTerm = coeff(relation, 0)
    candidateSet = Set(candidates)

    for bits in (128, 256, 512, 1024, 2048, 4096)
        B = ArbField(bits)
        xballs = B.(xs)
        yballs = B.(ys)
        tballs = B.(ts)
        pairwiseDisjoint(xballs) || continue
        pairwiseDisjoint(yballs) || continue
        pairwiseDisjoint(tballs) || continue

        matched = Tuple{Int, Int}[]
        certified = true
        for tball in tballs
            if any(
                guard -> contains_zero(evaluate(guard, tball)),
                projection.guards,
            )
                certified = false
                break
            end

            denominator = evaluate(coefficientOfX, tball)
            if contains_zero(denominator)
                certified = false
                break
            end

            xball = -evaluate(constantTerm, tball) / denominator
            yball = (tball - xball) / B(lambda)
            xmatches = findall(ball -> overlaps(xball, ball), xballs)
            ymatches = findall(ball -> overlaps(yball, ball), yballs)
            if length(xmatches) != 1 || length(ymatches) != 1
                certified = false
                break
            end

            pair = (only(xmatches), only(ymatches))
            if !(pair in candidateSet) || pair in matched
                certified = false
                break
            end
            push!(matched, pair)
        end

        if certified && length(matched) == length(candidates)
            sort!(matched)
            return matched
        end
    end
    return nothing
end

# Solve H(x,y)=q(x)=0 by factoring q and matching only within each irreducible
# x-block.  Local y-resultants avoid the large global Cartesian product; a
# third exact projection certifies the remaining coordinate correspondence.
function pointsFromXProjection(H, q, precise::Int)
    q = univariateSquarefreePart(q)
    degree(q) <= 0 && return NamedTuple[]
    points = NamedTuple[]

    for (localFactor, _) in factor(q)
        xs = realRootsAsAlgebraic(localFactor)
        isempty(xs) && continue

        R = parent(H)
        localFactorXY = fromUnivariate(localFactor, R, 1)
        localYResultant = resultant(localFactorXY, H, 1)
        iszero(localYResultant) &&
            throw(ArgumentError("the system is not zero-dimensional"))
        localYPolynomial = toUnivariate(localYResultant, 2)
        localYPolynomial === nothing &&
            error("internal error: local resultant retained x")
        ys = realRootsAsAlgebraic(localYPolynomial)
        isempty(ys) && continue

        candidates = arbCandidatePairs(H, xs, ys)
        isempty(candidates) && continue

        matched = nothing
        for magnitude in 1:32
            for lambda in (magnitude, -magnitude)
                matched = certifyCandidatePairs(
                    H,
                    localFactor,
                    xs,
                    ys,
                    candidates,
                    lambda,
                )
                matched === nothing || break
            end
            matched === nothing || break
        end
        if matched === nothing
            # Degenerate projections are rare.  Retain the original exact
            # substitution as a local fallback, after factorization and Arb
            # rejection have already made the candidate set small.
            matched = [
                (i, j) for (i, j) in candidates
                if iszero(evaluateAtAlgebraicPoint(H, xs[i], ys[j]))
            ]
        end

        for (i, j) in matched
            a = xs[i]
            b = ys[j]
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

# Generic fallback for a zero-dimensional system which has no supplied
# univariate x-equation.  This retains the original exact matcher.
function genericZeroDimensionalPoints(F, G, precise::Int)
    tx = eliminationPolynomial(F, G, 1)
    ty = eliminationPolynomial(F, G, 2)
    xs = realRootsAsAlgebraic(tx)
    ys = realRootsAsAlgebraic(ty)
    points = NamedTuple[]

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

# Return all distinct exact real solutions of F = G = 0.
function zeroDimensionalPoints(F, G, precise::Int)
    parent(F) == parent(G) ||
        throw(ArgumentError("F and G must have the same parent"))
    iszero(F) && throw(ArgumentError("F must not be zero"))
    iszero(G) && throw(ArgumentError("G must not be zero"))

    fInX = toUnivariate(F, 1)
    gInX = toUnivariate(G, 1)
    if fInX !== nothing && degree(fInX) > 0 && gInX === nothing
        points = pointsFromXProjection(G, fInX, precise)
        points === nothing ||
            return points
    elseif gInX !== nothing && degree(gInX) > 0 && fInX === nothing
        points = pointsFromXProjection(F, gInX, precise)
        points === nothing ||
            return points
    end

    return genericZeroDimensionalPoints(F, G, precise)
end

samePoint(p, q) = p.x == q.x && p.y == q.y


# Curve-isolation algorithm

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
