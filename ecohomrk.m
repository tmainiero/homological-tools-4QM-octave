function H=ecohomrk(psi,N,varargin)

%Computes the rank/dimension of the entanglement cohomology groups of
%the state psi of an N-partite system.  The output is a dimension vector
%[H^0,H^1,...,H^{N-1}]

%varargin is one of the of the following

%1) varargin=dim , i.e. the optional dimension vector specifying the
%Hilbert space dimension of each primitive subsystem


if nargin==3
  dim=varargin{:};
else
  dim=repmat(2,1,N);
end

[End_A,proj_A,simplex]=comfunctor(psi,N,dim);

H=zeros(1,N);

%[Z,Im_1,Op_basis_1]=kerim(psi,N,0,dim,End_A,proj_A{1,:},simplex);

[Z,Im_1]=kerim(1,psi,N,0,dim);
H(1)=size(Z,2)-1;
%Dimension of H^0 = ker(coboundary) -1 = # of column vectors in the
%nullspace Z=ker(cobdry_0).  We subtract 1 because the image of C_{-1} =
%\mathbb{C} is the identity on each 0-simplex.

for n=1:N-1
  %[Z,Im_2,Op_basis_2]=kerim(psi,N,n,dim,End_A,simplex,Op_basis_1);
  [Z,Im_2]=kerim(1,psi,N,n,dim);
  H(n+1)=size(Z,2)-size(Im_1,2);
  %[size(Z,2),size(Im_1,2)]
  
  %Dimension of H^{n-1} is rank ker_{n-1} - rank Im_{n-2}.

  Im_1=Im_2;
end
