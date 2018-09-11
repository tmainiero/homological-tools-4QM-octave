function rho=wst(N,varargin)

%Computes the (unnormalized) "W-state" = |10....> + |01...> + \cdots for an N-partite system.  If varargin is null
%all systems are assumed to be qubit systems; otherwise varargin is
%the dimension vector for the collection of subsystems.

if nargin==2
  dim=varargin{:};
else
  dim = repmat(2,[1,N]);
end

dim_tot=prod(dim);
psi_mat=zeros([N,dim_tot]);
for k = 1:N
  psi_mat(k,:) = ket(eye(N)(k,:),dim);
end

psi=sum(psi_mat,1);

rho=psi'*psi;
