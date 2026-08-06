@testset "Exact isolation" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    circle = x^2 + y^2 - 1

    critical_points = getPointsCritical(circle)
    @test length(critical_points) == 2

    boxes = bivariateRealIsolation(circle, 16)
    @test length(boxes) == 2

    matched_boxes = criticalPointBoxes(
        circle;
        boxPrecision = 24,
        intervalPrecision = 12,
    )
    @test length(matched_boxes) == 2
end
