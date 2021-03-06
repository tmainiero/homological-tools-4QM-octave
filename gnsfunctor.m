function [End_A,proj_A,simplex]=gnsfunctor(psi,N,varargin)

% Computes the endomorphism group End(G_{A}) for each subset A, where
% G_{A}<H_{A} is the Hilbert space generated by non-zero eigenvectors of
% the subsystem state rho_{A}.  Also outputs the projection operators
% \rho_A \in End(H_{A}) which projects onto the subspace End(G_A) for
% sets of order 1.

max_sets=nchoosek(N,floor(N/2));
%Maximum number of subsets/simplices of a given order.

if nargin==3
  dim=varargin{:};
  max_dim=max(dim)^N;
  [rho_A,simplex]=substates(psi,N,dim);
  End_A=cell(N+1,max_sets,max_dim,max_dim);
else
 % dim=repmat(2,[1,N]);
  [rho_A,simplex]=substates(psi,N);
  max_dim=2^N;
  End_A=cell(N+1,max_sets,max_dim,max_dim);
end

proj_A=cell(N+1,max_sets);
%dim_red=cell(N+1);

for ord=0:N
  set_ind_max=nchoosek(N,ord);
  for set_ind=1:set_ind_max
    [V,D]=eig(rho_A{ord+1,set_ind});
    %Determine eigenvectors V, organized as column vectors, and diagonal
    %form D for rho_A
    
    L=find(abs(diag(D))>2*eps);
   %Find indices of all non-zero eigenvalues

   %dim_red{ord+1}(set_ind)=size(L,1);
    %Dimensions of G_{A}.  dim{ord+1} outputs a vector of dimensions
    %with the same lexicographical ordering of simplices.

   V=V(:,L);
   %Only keep eigenvectors corresponding to non-zero eigenvalues
   

   proj_A{ord+1,set_ind}=V*V';

   %Projection operator on subspace G_{A} is given by taking the
   %outer-product of all (norm 1) eigenvectors with themselves and summing.

   if size(L,1)>0
     for l=1:size(D,1)
       for m=1:size(L)
	 End_A{ord+1,set_ind,l,m}=eye(size(D,1))(:,l)*V(:,m)';
       end
     end
   else
     End_A{ord+1,set_ind,1,1}=0;
   end
  end
end
