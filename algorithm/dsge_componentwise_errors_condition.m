function [results] = dsge_componentwise_errors_condition(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%24/03/2022
A=inputs.A;
B=inputs.B;
C=inputs.C;
D=inputs.D;
P=inputs.P;
[ny,~]=size(P);
F=A*P+B;
P_F=norm(P,'fro');


R_P=A*P^2+B*P+C;
[temp_P_FE in difp]=SYLG_allocated(ny,ny,real(F),eye(ny),A,real(P)',real(R_P),1);
results(8,1)=difp;
results(14,1)=norm(temp_P_FE)/P_F;

gamma=@(k)k*eps/(1-k*eps);
    %tol=results(6,1)*ny^2*(eps+gamma(ny+2)+gamma(2*ny+2));
R_u=eps*abs(C)+gamma(ny+2)*abs(B)*abs(P)+gamma(2*ny+2)*abs(A)*abs(P)^2;
[temp_P_FE in2 difp2]=SYLG_allocated(ny,ny,real(F),eye(ny),A,real(P)',abs(R_P)+R_u,1);
results(8,2)=difp2;
results(14,2)=max(abs(temp_P_FE),[],"all")/max(abs(P),[],"all");