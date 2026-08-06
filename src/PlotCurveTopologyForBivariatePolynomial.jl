module PlotCurveTopologyForBivariatePolynomial

include("isolation/bivariate_real_isolation.jl")
include("isolation/qbar_critical_points.jl")
include("isolation/match_critical_boxes.jl")
include("topology/connect_critical_boxes.jl")
include("topology/count_plane_regions.jl")
include("plotting/plot_critical_box_edges.jl")

export bivariateRealIsolation
export getPointsCritical, getPointsCriticalWithIntervals
export matchBoxesToCriticalIntervals, criticalPointBoxes
export getTopologySampleBoxes, topologyGraphData, connectCriticalBoxes
export plotCriticalBoxEdges, plotCriticalBoxes
export compactifiedTopologyGraphData, countPlaneRegions

end
