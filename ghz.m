function rho=ghz(N,varargin)

%Creates the (unnormalized) GHZ (density) state for an N-partite qubit system (this
%is the Bell state |00> + |11> when N=2).  When varargin is present we
%generalize (by ``lifting": embed \mathbb{C}^2 into \mathbb{C}^N via
%the first two coordinates) this state to an N-partite system with
%dimension vector dim = varargin.

if nargin==2
  dim=varargin{:};
else
  dim=repmat(2,[1,N]);
end

psi=ket(repmat(0,[1,N]),dim) + ket(repmat(1,[1,N]),dim);

rho=psi'*psi;
