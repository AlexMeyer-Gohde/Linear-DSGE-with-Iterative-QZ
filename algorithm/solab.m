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

function [n,l,f,p,eig] = solab(a,b,phi,c,nk);

[s,t,q,z] = qz(a,b);                % upper triangular factorization of the matrix pencil b-za
eig=diag(t)./diag(s);
[s,t,q,z] = ordqz(s,t,q,z,'udo');   % reordering of generalized eigenvalues with the block inside the unit circle in the upper left
z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);
nd=length(z)-nk;
if rank(z11)<nk;
	error('Invertibility condition violated')
end

z11i = z11\eye(nk);
s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);

if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
   warning('Wrong number of stable eigenvalues.');
end

dyn = s11\t11;
f = real(z21*z11i);
p = real(z11*dyn*z11i);

z22=z(nk+1:end,nk+1:end);
z12=z(1:nk,nk+1:end);
t12=t(1:nk,nk+1:end);
s12=s(1:nk,nk+1:end);
nx=length(phi);


if nk<20 %With a "small" number of endogenous variables, direct calculation quicker
    TEMP_1=kron(phi',s(2*nk-nd+1:2*nk,2*nk-nd+1:2*nk))...
        -kron(eye(nx,nx),t(2*nk-nd+1:2*nk,2*nk-nd+1:2*nk));
    TEMP_2=q(2*nk-nd+1:length(q),:)*c;
    vec_m=TEMP_1\reshape(TEMP_2,numel(TEMP_2),1); %Equation 5.7 in Klein (2000)
    m=reshape(vec_m,nd,nx); %Undo column-vectorization
else %Otherwise recursive calculation (eqs 5.9-5.13 in Klein (2000)) quicker
    QC=q(2*nk-nd+1:length(q),:)*c;%define lower half of q*c
    for j=1:nk
        r_tran=QC(nk+1-j,:);%Starting at the end
        for i=nk+2-j:nk %and working forwards
            r_tran=r_tran+t(2*nk+1-j,i+nk)*m(i,:)-s(2*nk+1-j,i+nk)*m(i,:)*phi;
        end
        m(nk+1-j,1:nx)=r_tran/(s(2*nk+1-j,2*nk+1-j)*phi-t(2*nk+1-j,2*nk+1-j)*eye(nx));
    end
end

n=real((z22-z21*z11i*z12)*m);
l=real(-z11*dyn*z12*m+z11*s11\(t12*m-s12*m*phi+q(1:2*nk-nd,:)*c)+z12*m*phi);
