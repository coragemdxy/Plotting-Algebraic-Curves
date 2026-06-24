function main()
    R,(x,y) = polynomial_ring(QQ, ["x", "y"])
    f = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*x
    f1 = 0
    for (g,e) in zip(coefficients(f), exponent_vectors(f))
        if e[2] == 2
            f1 += g*x^(e[1])
            println(f1)
        end
    end
    println(f1)
end