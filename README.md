License: GPLv2 or later

#Quickstart: Octave Cohomological Toolbox for the Quantum Mechanic

1. gnscohomrk(psi,N,dims)

2. ecohomrk(psi,N,dims)

Computes the dimensions of the GNS (case 1) or commutant (case 2) cohomology vector spaces of the state psi (pure or mixed) of an N-partite system.  The output is a dimension vector of the form [H^0, H^1,..., H^{N-1}].
where H^{i} is the rank of the ith GNS/commutant complex.  This is equivalently the list of coefficients of the associated Poincare polynomial.

-'psi' can be either:
	--a row vector representing a pure state, or a density state: e.g. psi = [1,0,0,0]  
	--A density state (a self-adjoint square matrix), e.g. psi = [1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0,0].

-'N' is the number of primitive subsystems (the "partiteness").

-'dims' is an optional dimension vector specifying the Hilbert space dimension of each primitive subsystem.  Must be a vector of length N.  If this is left empty, it will be assumed all subsystems are dimension 2.

-The inputs 'psi', 'N', and 'dims' must be sensible, i.e. psi must be able to be represented as an honest N-partite (density) state on primitive subsystems with dimension vector [d_1,...d_N].  That is to be written as an 1 x 2^{N} vector or a 2^{N} x 2^{N} (square matrix) if 'dims' is empty, or a length d_1d_2...d_N row vector (square matrix) if dims=[d_{1}, ..., d_{N}].

ket(basis,dims)

Ket vector given the standard computational basis on a multipartite system: outputs the vector |basis(1)> \otimes |basis(2)> \otimes \cdots \otimes |basis(n)> = |basis(1) ... basis(n)>. As a row vector in octave.

-dims is an optional row vector specifying the Hilbert space dimension of each primitive subsystem. If left empty it is assumed that there are length(basis) subsystems of dimension 2. 

-basis is a list of non-negative integers [b_1,b_2,...b_N], where b_{i} runs from 0 to dims(i) - 1.  Here b_{k} represents the standard basis element |b_{i}> of \mathbb{C}^{dims(i)}  with a 1 in the (d_{k} - 1)^{th} position, and zeros elsewhere.


ghz(N,dims)

The (unnormalized) GHZ density state given by (|00....0> + |11...1>)(<00...0| + <11..1|) for an N-partite qubit system.  (This is an unnormalized density state associated to the Bell state |00> + |11> when N=2).

-'N' is the number of primitive subsystems.


-'dims' is an optional dimension vector of length N of the ambient system: assumed to be a list of 2's when empty.

tensor(a,b,...)

Returns the kronecker product of its arguments.  To produce the state psi = |01> for instance, we would write
psi = tensor([0,1],[1,0]).


Taken from Toby Cubitt.

#Acknowledgements:
Basic linear algebraic operations: partialtrace.m, tensor.m, and syspermute.m are taken from Toby Cubitt's Matlab software (c.f. http://www.dr-qubit.org/matlab.html) licensed under GPLv2.
