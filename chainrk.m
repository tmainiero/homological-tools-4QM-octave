function C=chainrk(psi,N,varargin)

%Outputs the chain group dimensions



if nargin==3
  dim=varargin{:};
else
  dim=repmat(2,1,N);
end

[End_A,proj_A,simplex]=sheaf(psi,N,dim);

C=zeros(1,N);

for n=0:N-1
  ord=n+1;
  set_ind_max=nchoosek(N,ord);
  C_partial=zeros(1,set_ind_max);
  for set_ind=1:set_ind_max
    C_partial(1,set_ind)=rank(proj_A{ord+1,set_ind})^2;
  end
  C(ord)=sum(C_partial,2);
end
