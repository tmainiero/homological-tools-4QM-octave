function [rho_A,simplex]=substates(psi,N,varargin)

% Takes in a total state psi (expressed as either a density state or a
% pure state) and the total number of subsystems N, and outputs the
% partial traces on all 2^N subsystems.  Third argument gives the
% dimensions of all subsystems.  Subsystems are ordered first by
% order, then by lexicographical ordering within each set of sets of
% fixed order (e.g. {1,2,3} < {1,2,4} < {1,3,4})

% The output is given as cell array rho(order,l) where 'order' indicates
% the set/simplex order and l indexes the lexicographical ordering of
% simplices of a given order.


if any(size(psi) == 1)
  % state vector
  if size(psi,1) == 1
    psi=psi';
  end
  rho=psi*psi';
  rho=rho/trace(rho);
  %normalize
else
  rho=psi;
  rho=rho/trace(rho);
  %normalize
end

if nargin==3
  dim=varargin{:};
else
  dim=repmat(2,1,N);
end
%Set dimensions of hilbert spaces of subsystems.  Default to N-qubit
%system if dimensions are not specified.

syslist=1:N;

rho_A{N+1,1}=rho;
%State after tracing over the empty set (the only order 0 set) is the
%empty set.

simplex=cell(N+1,1);

for co_ord=1:N
  ord=N-co_ord;
  tr_sys=nchoosek(syslist,co_ord);
  %List of subsystems to trace over; each row gives a subsystem of order
  %co_ord, which after tracing leaves a subsystem of order N-co_ord

  num_sets=size(tr_sys,1);
  set_ind=num_sets-(1:num_sets)+1;
  tr_sys=tr_sys(set_ind,:);

  %set_ind indexes the leftover systems in lexicographical order
  %(which happens to have reversed order from the traced systems), and
  %we permute the rows of tr_sys accordingly.


  for l=1:num_sets
    %l indexes the systems to be traced over/simplices in
    %lexicographical order

    simplex{ord+1}(l,:)=setdiff(syslist,tr_sys(l,:));
    %Outputs simplices corresponding to index l.  simplex{ord+1}(l,:)
    %provides the lth simplex of order ord.

    rho_A{ord+1,l}=partialtrace(rho,tr_sys(l,:),dim);
  end
  simplex{N+1}=1:N;
  %Repairs the case when tr_sys=[]; so l=0 and the loop is not run.
end