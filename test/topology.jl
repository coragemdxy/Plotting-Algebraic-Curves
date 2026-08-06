@testset "Topology graph" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])

    circle_graph = topologyGraphData(
        x^2 + y^2 - 1;
        boxPrecision = 24,
        intervalPrecision = 12,
        samplesPerStrip = 1,
    )
    @test length(circle_graph.vertices) == 4
    @test length(circle_graph.edges) == 4

    isolated_graph = topologyGraphData(
        x^2 + y^2;
        boxPrecision = 24,
        intervalPrecision = 12,
        samplesPerStrip = 1,
    )
    @test length(isolated_graph.vertices) == 1
    @test isempty(isolated_graph.edges)

    @test !isempty(methods(plotCriticalBoxes))
end
