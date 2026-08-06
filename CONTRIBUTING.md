# Contributing

Thank you for helping improve `PlotCurveTopologyForBivariatePolynomial.jl`.

## Development setup

1. Install Julia 1.10 or later.
2. Clone the repository and instantiate the project:

   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. Run the test suite:

   ```bash
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

4. Build the documentation when changing public behavior or documentation:

   ```bash
   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
   julia --project=docs docs/make.jl
   ```

## Pull requests

- Keep each pull request focused on one change.
- Add tests for bug fixes and new behavior.
- Preserve exact arithmetic in every topology decision.
- Explain changes to public interfaces in the pull request and changelog.
- Do not commit generated documentation, LaTeX build files, or local plots.

The existing camelCase public names are retained for compatibility. New
internal helpers should follow established Julia naming conventions and use
`!` for mutating functions.

## Reporting issues

Use the provided GitHub issue forms. Include a minimal polynomial, Julia
version, package status, precision arguments, and the complete error message.
