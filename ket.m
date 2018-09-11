function psi=ket(expr,varargin)

%Takes in an expression expr = [1,0,1,0,2] and a dimension vector
%varargin=dim and outputs the standard ket vector associated to this expression.

N=size(expr,2);
%number of systems, inferred from expr

if nargin==2
  dim=varargin{:};
else
  dim=repmat(2,[1,N]);
end


if ~isequal(N,size(dim,2))
   error('Size of input does not match the number of systems in the dimension vector');
end

ind=expr+1;
%Index used in loop.  I.e. 0 --> [1,0,0,0,...], 1 --> [0,1,0,0,0,...],
%2 ---> [0,0,1,0,0....], ..., k --> vector with 1 in (k+1)th position
%and zeros elsewhere.

psi=1;
%Initialization

for k = 1:N
  nextsys=zeros([1,dim(k)]);
  nextsys(ind(k))=1;
  psi=tensor(psi,nextsys);
end
