function delta_Op_A=ecobdry(Op_A,proj_A,n,N,dim,varargin);

%Takes in an n-cochain Op_A, which assigns an operator in the space of
%sections End_A for each n-simplex A (an n-simplex consists of
%(n+1)-sets), the projection operators proj_A{ord,k} for every
%|A_{k}|=ord,  which act as identity operators on each End_A< End H_A,
%and outputs the coboundary: an (n+1)-cochain. Op_A is input as an array
%of matrices (Op_A{set_ind} is a matrix) and it is assumed that the sets of order n are ordered lexicographically.

%dim is a vector of size N which provides the dimensions of all
%subspaces H_A for all 0-simplices. Note: this is redundant information
%given purely Op_A.

%The last argument is an optional argument providing pre-computed
%simplices of all orders.

if n>=N-1
  delta_Op_A=0;
  %The N-th chain group is defined as the trivial vector space 0.
else

  ord=n+1;

  if nargin==4
    syslist=1:N;
    co_ord_1=N-ord;
    tr_sys_1=nchoosek(syslist,co_ord_1);
    tr_sys_2=nchoosek(syslist,co_ord_1-1);
    
    for l=1:size(tr_sys_1,1)
      set_ind=size(tr_sys_1,1)-l+1;
      simplex_1(set_ind,:)=setdiff(syslist,tr_sys_1(l,:));
    end
    
    if size(tr_sys_2,1)>0
      
      for l=1:size(tr_sys_2,1)
	set_ind=size(tr_sys_2,1)-l+1;
	simplex_2(set_ind,:)=setdiff(syslist,tr_sys_2(l,:));
      end
      
    else
      simplex_2=syslist;
    end
 %Corrects for situation where the n=N-1, i.e. top simplices so tr_sys_2=[];

  else
    simplex_1=varargin{:}{ord+1,:};
    %List of n-simplices

    simplex_2=varargin{:}{ord+2,:};
    %List of (n+1)-simplices
  end

  for l=1:size(simplex_2,1)

    bdry=nchoosek(simplex_2(l,:),ord);
    %Find all n-simplices inside of the lth n+1-simplex.
  
    num_simp=size(bdry,1);
    set_ind=num_simp-(1:num_simp)+1;
    bdry=bdry(set_ind,:);
    %Re-order backwards to be consistent with the boundary operator acting
    %on each simplex (e.g. [1,2,3]--> [2,3;1,3;1,2]).
  
    sys_cotr=zeros(num_simp,1);
    for k=1:num_simp
      sys_cotr(k)=setdiff(simplex_2(l,:),bdry(k,:));
    end
   %Find sets eliminated in simplex_2(l,:) to produce bdry These will be
   %used in taking cotraces (tensoring with identity).

    dim_cotr=dim(sys_cotr);
    %Find dimensions of Hilbert spaces G_{A} in sys_cotr

    if n>0
      [tf,bdry_row]=ismember(bdry,simplex_1,'rows');
    else
      [tf,bdry_row]=ismember(bdry,simplex_1);
    end
    %Find the row-location of the simplices of bdry inside of simplex_1.
    %The first argument is a useless logical (it should be a column vector
    %of 1's if everything goes well).
  
    dim_simp_2=prod(dim(bdry(1,:)))*dim_cotr(1);
    partial_delta=zeros(dim_simp_2,dim_simp_2,num_simp);
    %initialize

    for k=1:num_simp
      DIM=[dim(bdry(k,:)),dim_cotr(k)];
      %vector of dimensions of subsystems for kth partial boundary
    
      perm=[1:k-1,n+2,k:n+1];

      %permutation to perform after tensoring with identity on the right:
      %after performing the permutation we should have tensored with the
      %appropriate projection/identity in the kth position.
    
      if Op_A{bdry_row(k)}==0
	partial_delta(:,:,k)=zeros(dim_simp_2);
      else
	partial_delta(:,:,k)=(-1)^(k-1)*syspermute(kron(Op_A{bdry_row(k)},proj_A{2,sys_cotr(k)}),perm,DIM);
      end

    end
    %Tensor with identity, add appropriate sign, permute systems to
    %place identity in the location of the set eliminated to get
    %bdry(k,:) from simplex_2(l,:), then sum all partial boundaries and
    %project onto appropriate subspace.

    delta_Op_A{l}=proj_A{ord+2,l}*sum(partial_delta,3)*proj_A{ord+2,l};
  end

end
