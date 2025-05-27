function [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic)
if ~isempty(matrix_quadratic.P)
matrix_quadratic.P=matrix_quadratic.P(:,matrix_quadratic.index_m);
if matrix_quadratic.nstatic > 0
    temp = - matrix_quadratic.A_static*matrix_quadratic.P(matrix_quadratic.npred+1:end,:)*matrix_quadratic.P(matrix_quadratic.index_m,:);
    temp(:,matrix_quadratic.index_m) = temp(:,matrix_quadratic.index_m)-matrix_quadratic.C_static;
    temp = matrix_quadratic.B_static\(temp-matrix_quadratic.B_rest*matrix_quadratic.P(1:end,:));
    matrix_quadratic.P = [temp; matrix_quadratic.P];
    temp = [];
end
%matrix_quadratic.P=matrix_quadratic.P(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
matrix_quadratic.P=[zeros(matrix_quadratic.endo_nbr,matrix_quadratic.nstatic) matrix_quadratic.P zeros(matrix_quadratic.endo_nbr,matrix_quadratic.nfwrd)];
%matrix_quadratic.P=matrix_quadratic.P(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
end