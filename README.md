# DipoleChain.jl

Exact sparse diagonalization and path integral matrix multiplication for a linear chain of evenly-spaced [linear rigid rotors](https://en.wikipedia.org/wiki/Rigid_rotor#Quantum_mechanical_linear_rigid_rotor) with [dipole-dipole interactions](https://en.wikipedia.org/wiki/Intermolecular_force#Dipole-dipole_interactions).

Tested with Julia 0.6.


## Installation

1. `Pkg.clone("https://github.com/0/DipoleChain.jl.git")`


### Requirements

These should be pulled in automatically when installing this package.
To use it without installing it (e.g. from a local git checkout), you'll need to manually obtain the following dependencies:

* ArgParse (`Pkg.add("ArgParse")`)


## Units

The physical Hamiltonian `\hat{H} = B \hat{K} + (\mu^2 / 4 \pi \epsilon_0 r^3) \hat{V}` has several parameters:

* `B`: rotational constant (energy)
* `\mu`: [dipole](https://en.wikipedia.org/wiki/Dipole) moment
* `r`: separation between adjacent rotors (length)

The rotational energy operator `\hat{K}` (total squared angular momentum divided by `\hbar^2`) and the potential energy operator `\hat{V}` (dipole-dipole interaction between all rotor pairs) are dimensionless.

To simplify the situation, we instead use the one-parameter dimensionless Hamiltonian `\hat{H} / B = \hat{K} + \hat{V} / R^3`.
The only parameter is the separation distance `R` in natural length units, which may be found from the physical parameters using `R = r / (\mu^2 / 4 \pi \epsilon_0 B)^{1/3}`.
The resulting eigenvalues are in natural energy units and may be converted to physical energies by multiplying in the rotational constant `B`.

Reciprocal temperatures are in units of reciprocal energy (i.e. the [Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant) is set to 1).


## Examples

* `julia examples/diagonalization.jl --help`
* `julia examples/diagonalization.jl -R 1.23 -N 4 --l-max 3 --A-start 2 --A-size 1`
* `julia examples/path_integral.jl --help`
* `julia examples/path_integral.jl -R 1.23 -N 2 --l-max 3 --beta 45.67 -P 8`
* `julia examples/path_integral.jl -R 1.23 -N 2 --l-max 3 --tau 4.54 -P 3 --pigs --A-start 2 --A-size 1`


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
