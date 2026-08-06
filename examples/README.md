# Examples

The example files contain the original interactive `main` functions moved out
of package source. Load a file and call its `main` function explicitly, for
example:

```bash
julia --project=. -e 'include("examples/topology/count_plane_regions.jl"); main()'
```

The large plotting example is computationally expensive and writes its output
next to the example file.
