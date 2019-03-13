# LiquidsStructure.jl

[![Build Status](https://app.codeship.com/projects/9fcab2f0-271b-0137-ca37-1e86d73d396b/status?branch=master)

This library intended to provide a mean to compute the structure factor of a
variety of liquids with different interaction potentials and under different
approximation schemes (e.g. Percus–Yevick closure for the Ornstein–Zernike
relation for a hard-sphere liquid) in the [Julia](http://julialang.org)
programming language.

## Status

For the time being this library provides routines to calculate the structure
factor of:

 - hard-spheres under the Percus–Yevick closure for the OZ relation,
 - hard-spheres under the Verlet–Weis approximation,
 - hard-disks under the Rosenfeld FMT approximation,
 - dipolar hard-spheres under the MSA approximation.

**TODO**

- [ ] Documentation and examples
- [ ] Add more interaction potentials and approximation schemes

## Installation

`LiquidsStructure.jl` should work on Julia 1.0 and later versions and can be
installed from a Julia session by running

```julia
julia> Pkg.add(PackageSpec(url = "https://github.com/LANIMFE/LiquidsStructure.jl.git"))
```

## Usage

Once installed, run

```julia
using LiquidsStructure
```

The structure factor can be calculated in two ways:

 1. By calling the function `structure_factor(::InteractionPotential,
    ::ApproximationScheme, k)`, where `k` is a wavevector.
 -  By constructing a `StructureFactor(::InteractionPotential,
    ::ApproximationScheme)` object `S`, and then using `S` a function over the
    wavevector `k` (`S(k)`).
