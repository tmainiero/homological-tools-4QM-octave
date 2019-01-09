# Octave Cohomological Toolbox for the Quantum Mechanic

## License: GPLv2 or later
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Description

This is a set of Octave functions built to compute the cohomology of the cochain complexes introduced in the paper 'Homological Tools for the Quantum Mechanic' ([arXiv:1901.02011](https://arxiv.org/abs/1901.02011)).

For software that can compute ranks and explicit bases of cohomology components in Mathematica, see: <https://github.com/tmainiero/homological-tools-4QM-mathematica>.

## How to Download

### Git

`git clone https://github.com/tmainiero/homological-tools-4QM-octave.git`

### File by File

#### From the Github web interface:
1. Go to the file you want to download and click it to view the contents
2. Locate the "Raw" button (On the top right at the time of writing) and right click.
3. Save as...


# Descriptions of Basic Functions

## To compute a list of ranks of cohomology groups
1. `gnscohomrk(psi,N,dims)`

2. `ecohomrk(psi,N,dims)`

Computes the dimensions of the GNS (case 1) or commutant (case 2) cohomology vector spaces of the state psi (pure or mixed) of an N-partite system.  The output is a dimension vector of the form [H^0, H^1,..., H^(N-1)].
where H^(i) is the rank of the ith GNS/commutant complex.  This is equivalently the list of coefficients of the associated Poincare polynomial.

  * `psi` can be either:
    1. A row vector representing a pure state, or a density state: e.g. `psi = [1,0,0,0]` 
	2. A density state (a self-adjoint square matrix), e.g. `psi = [1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0,0]`.

  * `N` is the number of tensor factors/primitive subsystems (the "partiteness").

  * `dims` is an optional argument giving the *dimension vector*: the list of dimensions of Hilbert spaces at the tensor factors. It must be a vector of length N.  If this argument is left empty, it will be assumed all subsystems are dimension 2.

The inputs `psi`, `N`, and `dims` must be sensible, i.e. psi must be able to be represented as an honest N-partite (density) state on primitive subsystems with dimension vector [d_1,...d_N]: i.e. if `dims` is empty: it must be written as either an 1 x 2^(N) (row) vector or a 2^(N) x 2^(N) square matrix; if `dims=[d_1, ..., d_N]`, then it must be a length d_1d_2...d_N row vector or a d_1d_2...d_N square matrix.

## Basic State Manipulation

`ket(basis,dims)`

Ket vector given the standard computational basis on a multipartite system: outputs the vector |basis(1)> \otimes |basis(2)> \otimes \cdots \otimes |basis(n)> = |basis(1) ... basis(N)>. Output is a row vector in Octave.

  * `dims` is an optional row vector specifying the Hilbert space dimension of each primitive subsystem. If this argument is left empty it is assumed that there are N (= length of `basis`) subsystems of dimension 2. 

  * `basis` is a list of non-negative integers [b_1,b_2,...b_N], where b_i runs from 0 to dims(i) - 1.  Here b_k represents the standard basis element |b_i> of \mathbb{C}^(dims(i))  with a 1 in the (d_k - 1)^(th) position, and zeros elsewhere.

---

`ghz(N,dims)`

The (unnormalized) GHZ density state given by (|00....0> + |11...1>)(<00...0| + <11..1|) for an N-partite qubit system.  (This is an unnormalized density state associated to the Bell state |00> + |11> when N=2).

-`N` is the number of primitive subsystems.


-`dims` is an optional dimension vector of length N of the ambient system: assumed to be a list of 2's when empty.

---

`wst(N,dims)`

The (unnormalized) W density state given by (|0....01> + |0...10> + ...+ |1...0>)(<|0....01| + <0...10| + ...+ <1...0|) for an N-partite qubit system.  (This is an unnormalized density state associated to the Bell state |01> + |10> when N=2).

-`N` is the number of primitive subsystems.


-`dims` is an optional dimension vector of length N of the ambient system: assumed to be a list of 2's when empty.

---

`tensor(a,b,...)`

Returns the kronecker product of its arguments.  To produce the state psi = |01> for instance, we would write
`psi = tensor([0,1],[1,0])`.

This function is taken from Toby Cubitt.

# Acknowledgements:
Basic linear algebraic operations: partialtrace.m, tensor.m, and syspermute.m are taken from Toby Cubitt's Matlab software (c.f. http://www.dr-qubit.org/matlab.html) licensed under GPLv2.
