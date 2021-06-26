# LiquidsStructure.jl

[![GitHub Actions](https://github.com/LANIMFE/LiquidsStructure.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/LANIMFE/LiquidsStructure.jl/actions?query=workflow%3ACI)

This library intended to provide a mean to compute the structure factor of a
variety of model liquids with different interaction potentials and under
different approximation schemes (e.g. Percus–Yevick closure for the
Ornstein–Zernike relation for a hard-sphere liquid) in the
[Julia](http://julialang.org) programming language.

## Status

For the time being this library provides routines to calculate the structure
factor of:

 - hard spheres under the [Percus–Yevick
   closure](https://en.wikipedia.org/wiki/Percus–Yevick_approximation) for the
   [Ornstein–Zernike (OZ)
   equation](https://en.wikipedia.org/wiki/Ornstein–Zernike_equation),
 - hard spheres under the Percus–Yevick closure with [Verlet–Weis
   corrections](https://doi.org/10.1103/PhysRevA.5.939),
 - hard disks under the Rosenfeld FMT approximation,
 - dipolar hard spheres under the MSA approximation.

**TODO**

- [ ] Documentation
- [ ] Add more interaction potentials and approximation schemes

## Installation

`LiquidsStructure.jl` is compatible with Julia 1.0 and later versions. It requires first
adding the [LANIMFE-Registy](https://github.com/LANIMFE/LANIMFE-Registry) to your Julia
installation. Then it can simply be installed by running

```julia
julia> ]
pkg> add LiquidsStructure
```

## Usage

Once installed, run

```julia
using LiquidsStructure
```

The structure factor can be calculated by constructing a
`StructureFactor(::Liquid, ::ApproximationScheme)` object `S`, and then using
`S` as a function over the wavevector `k` (`S(k)`).

As an example, let us plot the structure factor for a hard spheres liquid with
a volume fraction `η = 0.4`, under the Verlet–Weis approximation scheme

```julia
julia> using LiquidsStructure
julia> using Plots
julia> η = 0.4;
julia> S = StructureFactor(HardSpheres(η), VerletWeis);
julia> plot(k -> S(k), 0, 10π, label = "S(k)")
```

![Example image](assets/example.png?raw=true)

## Acknowledgements

This project was developed with support from CONACYT through the Laboratorio
Nacional de Ingeniería de la Materia Fuera de Equilibrio (LANIMFE).
