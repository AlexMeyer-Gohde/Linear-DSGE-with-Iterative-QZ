function [diff] = convergence_criterion_dynare(M,P,P0,Bplus,A,Bminus,C,method)
% diff returns the convergence criterion
%
% INPUTS
%   M, Pj, P0, G,
%   A, B, C      [matrices]     matrices M,Pj,P0,G,A,B,C                               
%
%   method       [string]       A string determining which method is used 
%                               for calculating the convergence criterion   
%
% OUTPUTS
%   diff         [scalar]       The convergence criterion
%
%
% Copyright (C) 2022 Alexander Meyer-Gohde & Johanna Saecker
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% the default convergence criterion follows ...
if strcmp(method,'residual')
    diff=norm(M,'fro');%=norm((A*P+B)*P+C,'fro');
elseif strcmp(method,'diff')
    diff=norm(P-P0,1);
elseif strcmp(method,'reldiff')
    diff=norm(P-P0,1)/norm(P,1);
% elseif strcmp(method,'reldiffX')
%     diff=norm(P-P0,1)/norm(X,1);
elseif strcmp(method,'fe1')
    [ndynamic,nspred]=size(P);
    [ndynamic,nsfwrd]=size(Bplus);
    [ndynamic,nspred]=size(Bminus);
    nboth=nspred+nsfwrd-ndynamic;
    npred=nspred-nboth;
    nfwrd=nsfwrd-nboth;
    P_F=norm(P,'fro');
    Pplus=P(npred+1:end,:);
    Pminus=P(1:nspred,:);
    R_P=(A*Pplus+Bminus)*Pminus+Bplus*Pplus+C;
    %F=A*P+B;
    % [temp_P_FE in difp]=SYLG_allocated(ny,ny,real(F),eye(ny),A,real(P)',real(R_P),1);
    % diff1=norm(temp_P_FE)/P_F;
    gamma=@(k)k*eps/(1-k*eps);
    Ru=eps*abs(C)+gamma(ndynamic+2)*abs(Bplus)*abs(Pplus)+gamma(ndynamic+2)*abs(Bminus)*abs(Pminus)+gamma(2*ndynamic+2)*abs(A)*abs(Pplus)*abs(Pminus);
    [temp_P_FE in difp]=SYLG_allocated(ny,ny,real(F),eye(ny),A,real(P)',abs(real(R_P))+Ru,1);
    diff2=max(abs(temp_P_FE),[],"all")/max(abs(P),[],"all");
    %[diff1 diff2]
    diff=diff2;
elseif strcmp(method,'fe1_sparse')
% R_P=M;
% V=kron(speye(size(B,1)),G)+kron(P',A);
%     temp_P=lsqminnorm(V,R_P(:));
% diff=norm(temp_P)/P_F;
else
    disp('Unknown Convergence Criterion');
end
end

