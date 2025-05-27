%
% Function: solab
%
% Written by Paul Klein
%
% Rewritten in November 2007 after a suggestion by Alexander Meyer-Gohde
%
% Format: [f,p] = solab(a,b,nk);
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nk
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: qz, ordqz
%
% Modified for ALGORITHMS
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

function [Delta_P,eig,varargout] = solab_P_dyn(a,b,nk);

qz_zero_threshold=1e-06;
qz_criterium=1-qz_zero_threshold;

[s,t,q_star,z] = qz(a,b, 'real');                % upper triangular factorization of the matrix pencil b-za
eig = ordeig(s, t);
select = abs(eig) > qz_criterium;
    sdim = sum(select);
    [s, t, ~, z] = ordqz(s, t, q_star, z, select);
    eig = ordeig(s, t);

% [s,t,q_star,z] = qz(a,b);                % upper triangular factorization of the matrix pencil b-za
% eig=diag(t)./diag(s);
% [s,t,q_star,z] = ordqz(s,t,q_star,z,'udo');   % reordering of generalized eigenvalues with the block inside the unit circle in the upper left
z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);
nd=length(z)-nk;
if rank(z11)<nk;
	error('Invertibility condition violated')
end

% options_.qz_zero_threshold=1-06;
% options_.qz_criterium=1+options_.qz_zero_threshold;
% [ss, tt, w, sdim, dr.eigval, info1] = mjdgges(a, b, options_.qz_criterium, options_.qz_zero_threshold);


if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
   warning('Wrong number of stable eigenvalues.');
end

Delta_P = real(z21/z11);

if nargout>3
    disp('Too many output arguments')
    return;
elseif nargout ==3
    output.s11=s(1:nk,1:nk);
    output.t11=t(1:nk,1:nk);
    output.z11=z(1:nk,1:nk);
    varargout{1}=output;
end
