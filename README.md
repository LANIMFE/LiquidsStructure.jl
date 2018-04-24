# LiquidsStructure.jl

This library intended to provide a mean to compute the structure factor of a
variety of liquids with different interaction potentials and under different
approximation schemes (e.g. Percus--Yevick closure for the Ornstein--Zernike
relation for a hard-sphere liquid) in [Julia](http://julialang.org)

## Status

For the time being this library provides routines to calculate the structure
factor of:

 - hard-spheres under the Percus--Yevick closure for the OZ relation,
 - hard-spheres under the Verlet--Weis approximation,
 - hard-disks under the Rosenfeld FMT approximation,
 - dipolar hard-spheres under the MSA approximation.

**TODO**

- [ ] Add more interaction potentials and approximation schemes

## Installation

`LiquidsStructure.jl` should work on Julia 0.6 and later versions and can be
installed from a Julia session by running

```julia
julia> Pkg.clone("https://gitlab.com/NE-SCGLE/LiquidsStructure.jl.git")
```

## Usage

Once installed, run

```julia
using LiquidsStructure
```

The structure factor can be calculated in two ways:

 - 1. By using the function `structure_factor(::InteractionPotential,
      ::ApproximationScheme, k)`, where `k` is a wavevector.

 - 2. By constructing a `StructureFactor(::InteractionPotential,
      ::ApproximationScheme)` object `S`, and then using `S` a function over
      the wavevector `k` (`S(k)`).
