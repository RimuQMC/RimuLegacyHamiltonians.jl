# RimuLegacyHamiltonians

This repository contains legacy Hamiltonians for Rimu.jl. These Hamiltonians are already or will eventually be superseeded by more general Hamiltonian types.

### Installing RimuLegacyHamiltonians

The package is currently not registered. It can be installed from the repository on GitHub.com. It is meant to be used with Rimu.jl.
Hit the `]` key at the Julia REPL to get into `Pkg` mode and type
```
pkg> add "https://github.com/RimuQMC/RimuLegacyHamiltonians.jl" Rimu
```

### Usage

The package is now installed and can be imported with
```julia-repl
julia> using Rimu, RimuLegacyHamiltonians
```

### Exported types

```julia
BoseFS2C # address type with two bosonic components 
BoseHubbardMom1D2C # Hubbard chain in momentum space for two-component bosons
BoseHubbardReal1D2C # Hubbard chain in real space for two-component bosons
```

### Example

```julia-repl
julia> using RimuLegacyHamiltonians; import Rimu

julia> fs = BoseFS2C((1,2,1),(0,1,0))
BoseFS2C(BoseFS{4,3}(1, 2, 1), BoseFS{1,3}(0, 1, 0))

julia> h = BoseHubbardMom1D2C(fs)
BoseHubbardMom1D2C(BoseFS2C(BoseFS{4,3}(1, 2, 1), BoseFS{1,3}(0, 1, 0)); ua=1.0, ub=1.0, ta=1.0, tb=1.0, v=1.0)

julia> ep = Rimu.ExactDiagonalizationProblem(h)
ExactDiagonalizationProblem(
  BoseHubbardMom1D2C(BoseFS2C(BoseFS{4,3}(1, 2, 1), BoseFS{1,3}(0, 1, 0)); ua=1.0, ub=1.0, ta=1.0, tb=1.0, v=1.0),
  nothing;
  NamedTuple()...
)

julia> Rimu.solve(ep)
EDResult for algorithm LinearAlgebraSolver() with 15 eigenvalue(s),
  values = [-6.9403, -0.56511, -0.215232, 0.715199, 2.51409, 2.7328, 4.1605, 4.17157, 5.83875, 6.63714, 7.15727, 8.19517, 8.86239, 9.82843, 11.9073],
  and vectors of length 15.
  Convergence info: "Dense matrix eigensolver solution from `LinearAlgebra.eigen`", with howmany = 15 eigenvalues requested.
  success = true.

```