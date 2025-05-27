dynare_location='C:\dynare\6.2\matlab';
 save('dynare_location','dynare_location')

addpath('..\..\algorithm\')

clear all
 p = path; 
 load('dynare_location')
 addpath(dynare_location)
 dynare jermann98_moments_comparison noclearall;
 
 
 %%%%%%%%%%%%%%%%%%
gamma=@(k)k*eps/(1-k*eps);
ny=M_.endo_nbr;
options.convergence_tolerance=errors_dynare(10,1)*ny^2*(eps+gamma(ny+2)+gamma(2*ny+2));%errors_dynare(6,1)*ny^2*(eps+gamma(ny+2)+gamma(2*ny+2));
options.convergence_metric='fe1';
[P,output] = QZ_iterate(matrix_quadratic,options);
output.j
oo_improve.dr=oo_.dr;
oo_improve.dr.ghu=output.Q(oo_.dr.order_var,:);
oo_improve.dr.ghx=P(oo_.dr.order_var,oo_.dr.order_var);
oo_improve.dr.ghx=oo_improve.dr.ghx(:,nstatic+1:end-nfwrd);
option=options_;
option.varlist=M_.endo_names;
option.qz_criterium=1+eps;

[ivar,vartan,options_] = get_variables_list(option,M_);
moments_improve= th_autocovariances(oo_improve.dr,ivar,M_,option,option.nodecomposition);
%moments_improve=compute_model_moments(oo_improve.dr,M_,option);
oo_improve.var=moments_improve{1};

std_Y_improve=(oo_improve.var(Y_gr_index,Y_gr_index))^(1/2);
std_C_improve=(oo_improve.var(C_gr_index,C_gr_index)/oo_improve.var(Y_gr_index,Y_gr_index))^(1/2);
std_I_improve=(oo_improve.var(I_gr_index,I_gr_index)/oo_improve.var(Y_gr_index,Y_gr_index))^(1/2);
rp_improve=-400*oo_improve.var(R_index,M_index);

rf_improve=-400*(oo_.steady_state(M_index)+0.5*oo_improve.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_improve.dr.ghu(oo_.dr.inv_order_var(M_index),:)');

actual_improve=[rp_improve rf_improve std_Y_improve std_C_improve std_I_improve  ];

matrix_quadratic.X=P;
matrix_quadratic.P=P;
matrix_quadratic.Q=output.Q;
%quad_error=matrix_quadratic_backward_errors_short(matrix_quadratic);
[errors_improve] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_improve=matrix_quadratic;
 save('improve_results','targets','actual_improve','matrix_quadratic_improve','errors_improve')
