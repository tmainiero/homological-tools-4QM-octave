function [Z,Im,Op_basis_2]=kerim(func_type,psi,N,n,varargin)

%Computes the kernel and image of the coboundary map \delta_{n}: C_{n} -->
%C_{n+1}.  This function can be run on its own, in which case varargin
%is at most 1 extra argument, or can be called upon in the function
%ecohom, in which case varargin consists of several arguments so that
%certain quantities do not have to be recomputed in the computation of cohomology.

%Here varargin is one of the following:

%1) varargin=dim , i.e. the optional dimension vector specifying the
%Hilbert space dimension of each primitive subsystem

%2) varargin=[dim,End_A,proj_A,simplex]

%2) varargin=[dim,End_A,proj_A,simplex,Op_basis_1], here Op_basis_1 gives a
%basis for the nth chain group.

if func_type==0
  functor = @gnsfunctor;
  cobdry = @gnscobdry;
elseif func_type==1
  functor = @comfunctor;
  cobdry = @ecobdry;
else
    error('The first argument must either be 0---for the GNS functor---or 1, for the commutant functor');
end 

if nargin==5
  dim=varargin{:};
  [End_A,proj_A,simplex]=functor(psi,N,dim);
elseif nargin==4
  dim=repmat(2,1,N);
  [End_A,proj_A,simplex]=functor(psi,N);
elseif nargin==8
  dim=varargin{1};
  End_A=varargin{2};
  proj_A=varargin{3};
  simplex=varargin{4};
else
  if nargin ~=9
    error("Inappropriate number of arguments.")
  end
  dim=varargin{1};
  End_A=varargin{2};
  proj_A=varargin{3};
  simplex=varargin{4};
  Op_basis_1=varargin{5};
end
%Compute functor/cochain groups.

ord=n+1;
num_sets=nchoosek(N,ord);
%max_hilb_dim=max(dim)^(ord);

if n>N-1
  error("Trivially vanishing cohomology group.")

elseif n==N-1

  if nargin==8
    chain_1_dim=size(Op_basis_1,1);
  else
    num_nsimp_1=nchoosek(N,ord);
    
    End_reshape=reshape(End_A(ord+1,:,:,:),[size(End_A,2),size(End_A,3),size(End_A,4)]);
    ind_ne=find(~cellfun(@isempty,End_reshape));
    S=size(End_reshape);
    [set_ne,l_ne,m_ne]=ind2sub(S,ind_ne);
    %Finds non-empty indices of End_A (on sets |A|=ord) in order to run
    %loop efficiently.  ind_ne is the linear index of non-empty indices
    %and ind2sub converts these into actual subscripts.

    chain_1_dim=size(set_ne,1);
    %Number of non-empty indices = dimension of chain group C_{n}.

    Op_basis_1=cell(chain_1_dim,num_nsimp_1);
    %Initialize
  
    for b_ind=1:chain_1_dim
      k=set_ne(b_ind);
      l=l_ne(b_ind);
      m=m_ne(b_ind);
      Op_basis_1(b_ind,k)=End_A{ord+1,k,l,m};
    end
    empty_ind=cellfun(@isempty,Op_basis_1);
    Op_basis_1(empty_ind)=0;
 
  end
%Generates a basis for the chain group C_{N-1}=\oplus_{k}C_{N-1,k}: for
%each k, we select an operator which is a non-zero element on C_{N-1,k}
%and zero on all other summands.

  Op_basis_2=[];

  Z=eye(chain_1_dim);
  Im=[];


%The kernel (which is everything in the case where n=N-1), expressed as
%a matrix of vectors in the basis defined above.


else 
  if nargin==8
    chain_1_dim=size(Op_basis_1,1);
    
    num_nsimp_2=nchoosek(N,ord+1);

   End_reshape=reshape(End_A(ord+2,:,:,:),[size(End_A,2),size(End_A,3),size(End_A,4)]);
    ind_ne=find(~cellfun(@isempty,End_reshape));
    S=size(End_reshape);
    [set_ne,l_ne,m_ne]=ind2sub(S,ind_ne);
    %Finds non-empty indices of End_A (on sets |A|=ord) in order to run
    %loop efficiently.  ind_ne is the linear index of non-empty indices
    %and ind2sub converts these into actual subscripts.

    chain_2_dim=size(set_ne,1);
    %Number of non-empty indices = dimension of chain group C_{n}.

    Op_basis_2=cell(chain_2_dim,num_nsimp_2);
    %Initialize
  
    for b_ind=1:chain_2_dim
      k=set_ne(b_ind);
      l=l_ne(b_ind);
      m=m_ne(b_ind);
      Op_basis_2(b_ind,k)=End_A{ord+2,k,l,m};
    end
    empty_ind=cellfun(@isempty,Op_basis_2);
    Op_basis_2(empty_ind)=0;

%Generates a basis for the chain group C_{n+1}

  else
    num_nsimp_1=nchoosek(N,ord);  
    num_nsimp_2=nchoosek(N,ord+1);
    
    End_reshape=reshape(End_A(ord+1,:,:,:),[size(End_A,2),size(End_A,3),size(End_A,4)]);
    ind_ne=find(~cellfun(@isempty,End_reshape));
    S=size(End_reshape);
    [set_ne,l_ne,m_ne]=ind2sub(S,ind_ne);
    %Finds non-empty indices of End_A (on sets |A|=ord) in order to run
    %loop efficiently.  ind_ne is the linear index of non-empty indices
    %and ind2sub converts these into actual subscripts.

    chain_1_dim=size(set_ne,1);
    %Number of non-empty indices = dimension of chain group C_{n}.

    Op_basis_1=cell(chain_1_dim,num_nsimp_1);
    %Initialize
  
    for b_ind=1:chain_1_dim
      k=set_ne(b_ind);
      l=l_ne(b_ind);
      m=m_ne(b_ind);
      Op_basis_1(b_ind,k)=End_A{ord+1,k,l,m};
    end
    empty_ind=cellfun(@isempty,Op_basis_1);
    Op_basis_1(empty_ind)=0;

%Generates a basis for the chain group C_{n}=\oplus_{k}C_{n,k}: for each
%k, we select an operator which is a non-zero element on C_{n,k} and
%zero on all other summands.

  end

     End_reshape=reshape(End_A(ord+2,:,:,:),[size(End_A,2),size(End_A,3),size(End_A,4)]);
    ind_ne=find(~cellfun(@isempty,End_reshape));
    S=size(End_reshape);
    [set_ne,l_ne,m_ne]=ind2sub(S,ind_ne);
    %Finds non-empty indices of End_A (on sets |A|=ord) in order to run
    %loop efficiently.  ind_ne is the linear index of non-empty indices
    %and ind2sub converts these into actual subscripts.

    chain_2_dim=size(set_ne,1);
    %Number of non-empty indices = dimension of chain group C_{n}.

    Op_basis_2=cell(chain_2_dim,num_nsimp_2);
    %Initialize
  
    for b_ind=1:chain_2_dim
      k=set_ne(b_ind);
      l=l_ne(b_ind);
      m=m_ne(b_ind);
      Op_basis_2(b_ind,k)=End_A{ord+2,k,l,m};
    end
    empty_ind=cellfun(@isempty,Op_basis_2);
    Op_basis_2(empty_ind)=0;

%Generates a basis for the chain group C_{n+1}

  delta_Op_basis_1=cell(chain_1_dim,num_nsimp_2);
  for b=1:chain_1_dim
    Op_A=Op_basis_1(b,:);
    delta_Op_basis_1(b,1:num_nsimp_2)=cobdry(Op_A,proj_A,n,N,dim,simplex);
  end
%Construct array of coboundaries of basis n-chains


  delta_n_mat_part=zeros(chain_2_dim,chain_1_dim,num_nsimp_2);
  for r=1:chain_2_dim
    for c=1:chain_1_dim
      for k=1:num_nsimp_2
	delta_n_mat_part(r,c,k)=trace(Op_basis_2{r,k}'*delta_Op_basis_1{c,k});
      end
    end
  end

  delta_n_mat=sum(delta_n_mat_part,3);

%Computes \delta_{n} as a matrix C_{n} --> C_{n+1} with respect to the
%chosen basis.

  Z=null(delta_n_mat);
%Matrix of column vectors in the kernel of \delta_{n}.

  R=rref(delta_n_mat);
%Reduced row echelon form of delta_n_mat.  Row vectors span the image.
  
  Im=R(any(R,2),:)';
%Eliminate zero rows and transpose; image of delta_n_mat is presented
%via a basis of column vectors.

end

%[n,chain_1_dim]
%size(Op_basis_1,1);
