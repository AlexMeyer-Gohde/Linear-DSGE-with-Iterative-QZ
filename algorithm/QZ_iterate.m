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
A=matrix_quadratic.A;
B=matrix_quadratic.B;
C=matrix_quadratic.C;

if nargout>2
    disp('Too many output arguments')
    return;
elseif nargout ==2
    D=matrix_quadratic.D;
end


ny=size(A,1);


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
    P0=zeros(ny,ny);
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
    max_it=1000;%Default max_it
end





diff=tol+1;
j=0; 





while diff>tol % As long as the iterations haven't converged...

B_k=[ -A*P0 -A;  eye(ny,ny) zeros(ny,ny)];
A_k=[C+B*P0 B;   P0 eye(ny,ny)];
[Delta_P,eig] = solab_P(B_k,A_k,ny);
P=P0+Delta_P;
    diff=convergence_criterion([],P,P0,[],A,B,C,convergence_metric);           % determine convergence criterion
    %diff=norm(X-X_0,1)/norm(X,1);
    j=j+1;                               % advance the counter
    P0=P;
    if j>max_it
        disp('Maximum iterations reached')
        break
    end
end


% determine optional output
if nargout>2
    disp('Too many output arguments')
    return;
elseif nargout ==2
    output.Q=-(A*P+B)\D;
    output.j=j;
    varargout{1}=output;
end




