using Documenter
using PlotCurveTopologyForBivariatePolynomial

makedocs(
    modules = [PlotCurveTopologyForBivariatePolynomial],
    sitename = "PlotCurveTopologyForBivariatePolynomial.jl",
    format = Documenter.HTML(edit_link = "main"),
    pages = [
        "Home" => "index.md",
        "Algorithms" => "algorithms.md",
        "Public API" => "api.md",
    ],
    checkdocs = :none,
)
