#=
This program is used to realize the algorithm of real isolation of a univariate function, 
returning the interval containing a root.
=#
using Nemo

#Make a polynomial squarefree of an univariate polynomial
function squareFree(P)
    if iszero(P)
        throw(ArgumentError("zero polynomial has no squarefree part"))
    end

    dP = derivative(P)
    G = gcd(P,dP)
    return divexact(P, G)
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

#realisation of real isolation.
function realIsolation(P,precise)
    if degree(P) <=0
        return []
    end

    p = squareFree(P)
    M = boundOfRoots(p)

    interval = [-M-1, M+1]

    res = realIsolationPart(p,interval,precise,0)
    return res
end

function main()
    R,x = polynomial_ring(QQ, "x")
    f =x*(x - 2)*(x + 3)
    precise = 64
    println(realIsolation(f,precise))
end
