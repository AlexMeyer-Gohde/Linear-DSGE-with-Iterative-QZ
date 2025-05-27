function [P,varargout] = QZ_iterate(matrix_quadratic,varargin)
% Returns a (the minimal) solvent X of the matrix quadratic equation
% A P^2 +B P +C=0
%
% INPUTS
%   matrix_quadratic    [structure] A structure containing the square
%                                   matrices A, B, and C
%
%   options (optional)  [structure] A structure containing options
%
% OUTPUTS
%   P                   [matrix]    A solvent of 0=A*P^2+B*P+C
%
%   output (optional)   [structure] A structure containing the following
%                                   possible outputs:
%
%   diff (optional)     [scalar]    The last value of the convergence
%                                   criteron
%
%   j    (optional)     [scalar]    The number of iterations
%
%   resid (optional)    [matrix]    The residual of A*X^2+B*X+C
%
% ALGORITHMS
%   Johannes Huber and Alexander Meyer-Gohde (2025). 
%   QZ Iterate
%
%
% Copyright (C) 2025 Johannes Huber and Alexander Meyer-Gohde
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Initialization. Housekeeping: Extract the matrices from the struct
if matrix_quadratic.nspred>0
A=matrix_quadratic.AA;
B=matrix_quadratic.BB;
C=matrix_quadratic.CC;

if nargout>2
    disp('Too many output arguments')
    return;
elseif nargout ==2
    %D=matrix_quadratic.D;
end


nstatic=matrix_quadratic.nstatic;
nfwrd=matrix_quadratic.nfwrd;
npred=matrix_quadratic.npred;
nboth=matrix_quadratic.nboth;
nsfwrd=matrix_quadratic.nsfwrd;
nspred=matrix_quadratic.nspred;
ndynamic=matrix_quadratic.ndynamic;
index_m=matrix_quadratic.index_m;
npred=matrix_quadratic.npred;


% Check the number of input arguments
if nargin>2
    disp('Too many input arguments')
    return;
else    % allocate the correct options depending on the chosen algorithm
    if nargin==2
        options=varargin{1};
    end
end

% if nargin==2 && isfield(options,"method")
%     if strcmp(options.method,'SF1')
%         method=1;
%     elseif strcmp(options.method,'SF2')
%         method=2;
%     else
%         disp('Unrecognized method, continuing using SF1');
%         method=1;
%     end
% else
%     disp('No method selected, using SF1');
%     method=1;
%end

if nargin==2 && isfield(options,"P0")
    P0=options.P0;
else
    P0=zeros(nsfwrd,nspred);
end

if nargin==2 && isfield(options,"convergence_tolerance")
    tol=options.convergence_tolerance;
else
    tol=ny*eps;%Default convergence criterion
    %gamma=@(k)k*eps/(1-k*eps);
    %tol=results(6,1)*ny^2*(eps+gamma(ny+2)+gamma(2*ny+2)); %results(6,1)
    %is the condition number
end

if nargin==2 && isfield(options,"convergence_metric")
    convergence_metric=options.convergence_metric;
else
    convergence_metric="reldiff";%residual";
end

if nargin==2 && isfield(options,"max_it")
    max_it=options.max_it;
else
    max_it=100;%Default max_it
end





diff=tol+1;
j=0; 



A_zero_minus=[B(:,1:npred) 0.5*B(:,npred+1:nspred)];
A_zero_plus=[0.5*B(:,npred+1:nspred) B(:,nspred+1:ndynamic)];


eye_minus=[zeros(nboth,npred) eye(nboth)];
eye_plus=[eye(nboth) zeros(nboth,nfwrd)];

while diff>tol % As long as the iterations haven't converged...

B_k=[ -A_zero_minus-A*P0 -A;  eye_minus zeros(nboth,nsfwrd)];
A_k=[C+A_zero_plus*P0 A_zero_plus;   P0(1:nboth,:) eye_plus];
[Delta_P,eig,solout] = solab_P(B_k,A_k,nspred);
P=P0+Delta_P;
temp=solout.s11\solout.t11;
temp=temp/solout.z11;
temp=solout.z11(1:npred,:)*temp;
P_full=[real(temp);P];
    diff=convergence_criterion([],[P_full zeros(ndynamic,nfwrd)],P0,[],[zeros(ndynamic, npred) A],B,[C zeros(ndynamic,nfwrd)],convergence_metric);           % determine convergence criterion
    diff
    %diff=norm(X-X_0,1)/norm(X,1);
    j=j+1;                               % advance the counter
    P0=P;
    if j>max_it
        disp('Maximum iterations reached')
        break
    end
end


P=P_full;
else
    j=0;
    P=[];
end


% determine optional output
if nargout>2
    disp('Too many output arguments')
    return;
elseif nargout ==2
%    output.Q=-(A*P+B)\D;
    output.j=j;
    output.diff=diff;
    varargout{1}=output;
end




