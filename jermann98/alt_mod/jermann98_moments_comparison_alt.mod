
% Jermann's(1998) model by Hong Lan and Alexander Meyer-Gohde
%  - model is quarterly calibrated
%  - to annualize risky rate for example, we use  
%     - mean(r_yr) = 4*mean(r)
%     - std(r_yr)  = sqrt(4)*std(r)
%     - example of this annualization method: campell&cochrane (1999,p.218)
%  - CAC stands for capital adjustment cost

%--------------------------------------------------------------------------
% 0. Presets, please do not change
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% 2. Defining variables
%--------------------------------------------------------------------------
var k c i y  C_gr Y_gr I_gr // N N_gr    
    mu //m //r rf //rp erp     
    cost cost_der
    z 
      R M;//w d Rf rb Pb bp;

varexo e_z;

parameters
  tau b beta beta_star alpha a_bar
  rho_z 
  sigma_a_bar sigma_z_bar 
  cost_xi  cost_b  cost_c
  delta
  N_ss expi_expk
  Scac;

%--------------------------------------------------------------------------
% 3. Calibration
%--------------------------------------------------------------------------
//parameter values
load('Y_asset_test')
  tau       = X_0(1);%5;
  b         = X_0(2);%0.8175;%0.82;%0.8175;
  alpha     = 1 - 0.64;
  a_bar     = 0.005;
  beta_star = X_0(3);%0.993;%0.99;%0.9937;
  delta     = 0.025;
  cost_xi   = X_0(4);%4.4;%1/0.23;%4.3;
  rho_z       = 0.99;
  //sigma_a_bar = 0;
  sigma_z_bar =X_0(5);%0.009985;
  N_ss        = 1;   
 Scac = 1;     
//coefficient on adjustment cost function




%--------------------------------------------------------------------------
% 4. (Stationarized) Model 
%--------------------------------------------------------------------------
model;      
  //marginal utility of consumption
     mu = ( exp(c)-b*exp(c(-1))*exp(-(a_bar)) )^(-tau) 
           - b*beta*( exp(c(+1))*exp(a_bar)-b*exp(c) )^(-tau);
  //kernel
     exp(M) = beta*exp(-tau*(a_bar))*(mu/mu(-1));

  //capital accumulation
     exp(k) = (1-delta+Scac*cost )*exp(-a_bar)*exp(k(-1)) + (1-Scac)*exp(i);
  
  //production function
     exp(y) = exp(z)*exp(-alpha*(a_bar))*(exp(k(-1))^alpha)*(N_ss^(1-alpha));

  //labor equation
  //   1 - N_ss = 0;

  //market clearing
     exp(y) = exp(c) + exp(i);

  //capital adjustment cost: Jermann's (1998, p.273 ) specification 
    //adjustment cost function
       cost = (cost_b/(1-cost_xi))*( (exp(i)/exp(k(-1)))*exp(a_bar) )^(1-cost_xi) + cost_c;  
    //first derivative of adjustment cost function wrt investment/capital ratio
       cost_der = cost_b*( (exp(i)/exp(k(-1)))*exp(a_bar) )^(-cost_xi);    

  //wage equation
  //   exp(w) = (1-alpha)*exp(y)/N_ss; 

  //dividend
 //    exp(d) = exp(y) - exp(w)*N_ss - exp(i); 

  //risk free rate
  //   (1+(1/100)*rf)*(beta*exp(-tau*(a_bar+epsilon_a(+1)))*mu(+1)/mu) = 1;   

  //risky rate
exp(R) = ( alpha*(exp(a_bar))*(exp(y)/exp(k(-1))) 
                       + (Scac*cost+1-delta)/(Scac*cost_der+1-Scac)
                       - Scac*(exp(i)/exp(k(-1))*exp(a_bar)) )*(Scac*cost_der(-1)+1-Scac); 
// exp(R)=(1+(1/100)*r);
//exp(M)=m;%M=log(beta)-tau*(a_bar+epsilon_a)+log(mu)-log(mu(-1));%exp(M)=m;
//exp(Rf)=m(+1);%Rf=log(beta)-tau*(a_bar+epsilon_a(+1))+log(mu(+1))-log(mu);
  //bond price
   //  (1+(1/100)*rb) = (Pb+1)/Pb(-1); 

  //asset pricing equations
     //for risky rate
//        exp(R(+1))*beta*exp(-tau*(a_bar))*(mu(+1)/mu) = 1 ;
         exp(R(+1)+M(+1)) = 1 ;

     //for bond
   //     (1+(1/100)*rb(+1))*m(+1) = 1;

  //expected risk premium in percentage 
   //  erp = r(+1)-rf;
  
  //ex post risk premium in percentage
  //   rp = r - rf(-1);

  //ex post bond premium in percentage
  //   bp = rb - rf(-1);

  //growth rate in log, trend restored if any
     C_gr = c - c(-1) + a_bar;
     Y_gr = y - y(-1) + a_bar;
     I_gr = i - i(-1) + a_bar;
  //   N_gr = N - N(-1);
   
  //shock processes
    //stochastic growth trend
    //   epsilon_a = sigma_a_bar*e_a;
    //stochastic productivity level
       z         = rho_z*z(-1) + sigma_z_bar*e_z;
end;

%--------------------------------------------------------------------------
% 5. Calculating Steady State
%--------------------------------------------------------------------------
steady_state_model;

   expi_expk = 1-(1-delta)*exp(-a_bar);
   cost_c    = -cost_xi/(1-cost_xi)*(expi_expk*exp(a_bar));
   cost_b    = (expi_expk*exp(a_bar))^cost_xi;
  beta      = beta_star*(exp(a_bar))^(tau-1);
  cost      = (cost_b/(1-cost_xi))*( ( expi_expk )*exp(a_bar) )^(1-cost_xi) + cost_c;
  cost_der  = cost_b*( ( expi_expk )*exp(a_bar) )^(-cost_xi);
  m         = beta*exp(-tau*a_bar);
  k_N       = ( (1/m-(1-delta+Scac*cost)+expi_expk*exp(a_bar)*cost_der*Scac) / (alpha*exp((1-alpha)*a_bar)*(Scac*cost_der+1-Scac)) )^(1/(alpha-1));
  N         = log(N_ss);
  k         = log(k_N)+N;
  y         = alpha*k + (1-alpha)*N + (-alpha)*(a_bar);
  i         = log(expi_expk*exp(k));
  c         = log( exp(y)-exp(i) );
  mu        = ( exp(c)-b*exp(c)*exp(-a_bar) )^(-tau) - b*beta*( exp(c)*exp(a_bar)-b*exp(c) )^(-tau);
  r         = 100*( alpha*exp(a_bar)*exp(y)/exp(k)*(Scac*cost_der+1-Scac)+1-delta+Scac*cost-Scac*cost_der*expi_expk*exp(a_bar) -1);
  rf        = 100*(1/m-1);
  rb        = 100*(1/m-1);
  Pb        = 100/rb;
  epsilon_a = 0;
  z         = 0;  
  C_gr      = a_bar;
  Y_gr      = a_bar;
  I_gr      = a_bar;
  N_gr      = 0;
  rp        = 0;
  bp        = 0;
  erp       = 0;
  w         = log(1-alpha) + y - N;
  d         = log( exp(y)-exp(w)*exp(N)-exp(i)  ); 
R=log(1+(1/100)*r);
M=log(m);  
Rf=log(m);
end;

steady;

%--------------------------------------------------------------------------
% 6. Computation
%--------------------------------------------------------------------------
shocks;
 // var e_a = 1;
  var e_z = 1;
end;

%set_dynare_seed(81);
stoch_simul(noprint,periods=0, drop=0, irf=0, order=1);


C_gr_index=find(strcmp(M_.endo_names,'C_gr'));
Y_gr_index=find(strcmp(M_.endo_names,'Y_gr'));
I_gr_index=find(strcmp(M_.endo_names,'I_gr'));
M_index=find(strcmp(M_.endo_names,'M'));
R_index=find(strcmp(M_.endo_names,'R'));
Rf_index=find(strcmp(M_.endo_names,'Rf'));

std_Y=(oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_C=(oo_.var(C_gr_index,C_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
std_I=(oo_.var(I_gr_index,I_gr_index)/oo_.var(Y_gr_index,Y_gr_index))^(1/2);
rp=-400*oo_.var(R_index,M_index);

rf=-400*(oo_.steady_state(M_index)+0.5*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)*oo_.dr.ghu(oo_.dr.inv_order_var(M_index),:)');


100*(4*oo_.var(R_index,R_index))^(1/2);
100*(4*oo_.var(Rf_index,Rf_index))^(1/2);

targets=[6.18 0.8 0.01 0.51 2.65];
actual_dynare=[rp rf  std_Y std_C std_I];


create_matrix_quadratic_from_dynare
ALPHA_ZS_dynare=[zeros(n,nstatic) oo_.dr.ghx zeros(n,nfwrd)];
%ALPHA_ZS_dynare=ALPHA_ZS_dynare(oo_.dr.inv_order_var,oo_.dr.inv_order_var);%ALPHA_ZS_dynare=ALPHA_ZS_dynare(1:2,1:2);
matrix_quadratic.X=ALPHA_ZS_dynare;
BETA_ZS_dynare=oo_.dr.ghu; BETA_ZS_dynare=BETA_ZS_dynare(oo_.dr.inv_order_var,:);
matrix_quadratic.P=ALPHA_ZS_dynare;
matrix_quadratic.Q=BETA_ZS_dynare;
[errors_dynare] = dsge_backward_errors_condition_full(matrix_quadratic);
matrix_quadratic_dynare=matrix_quadratic;


save('dynare_results','targets','actual_dynare','matrix_quadratic_dynare','errors_dynare')

