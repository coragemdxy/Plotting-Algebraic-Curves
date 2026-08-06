@testset "Complement region counts" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])

    cases = [
        "circle" => (x^2 + y^2 - 1, 2),
        "two circles" => (
            (x^2 + y^2 - 1) * (x^2 + y^2 - 4),
            3,
        ),
        "line" => (y, 2),
        "hyperbola" => (x * y - 1, 3),
        "coordinate cross" => (x * y, 4),
        "isolated point" => (x^2 + y^2, 1),
    ]

    for (name, (polynomial, expected)) in cases
        @testset "$name" begin
            @test countPlaneRegions(
                polynomial;
                boxPrecision = 24,
                intervalPrecision = 12,
                samplesPerStrip = 1,
            ) == expected
        end
    end

    @test countPlaneRegions(zero(R)) == 0
end
