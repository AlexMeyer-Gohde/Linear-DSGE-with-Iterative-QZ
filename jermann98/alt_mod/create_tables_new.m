%%%%%
%macro
%3,8

%max P
%max eig

%rr
%BE
%mu

%eig sep
%sep
%psi
%phi

%fe1
%fe2
%%%%%
%macro short
%1,8

%max Q

%rr
%BE
%mu
%rr2
%BE2
%mu2

%sep
%psi1
%psi2
%psi3
%phi1
%phi2
%phi3

%fe1
%fe12
%fe2
%fe22

%%%%%
%macro
%3,8

%max PQ
%max eig

%rr
%BE
%mu

%eig sep
%sep
%psi
%phi

%fe1
%fe2
%%%%%
P_select=[2 4 5 9 7 10 12 14 15];
Q_select_row=[2 4 5 2 4 5 7 10 10 10 12 12 12 14 14 15 15];
Q_select_col=[2 2 2 3 3 3 2 2  3  4  2  3  4   3  4  3  4];
Q_select=sub2ind([15 5],Q_select_row,Q_select_col);
results_P=nan(15,8);
results_Q=nan(18,8);
results_PQ=nan(15,8);

load('aim_results'); cs=1; rp_exact=actual_aim(5);
results_P([1 3:6],cs)=targets([5 4 1 2 3])';
results_PQ([1 3:6],cs)=targets([5 4 1 2 3])';

load('solab_results');cs=2;  errors_w=errors_solab; actual_w=actual_solab;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('gensys_results');cs=3; errors_w=errors_gensys;actual_w=actual_gensys;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('uhlig_results');cs=4; errors_w=errors_uhlig;actual_w=actual_uhlig;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('dynare_results');cs=5; errors_w=errors_dynare;actual_w=actual_dynare;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('aim_results');cs=6; errors_w=errors_aim;actual_w=actual_aim;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('bp_results');cs=7; errors_w=errors_bp;actual_w=actual_bp;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

load('dynare_cr_results');cs=8; errors_w=errors_dynare_cr;actual_w=actual_dynare_cr;
results_P([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_P(2,cs)=rp_exact-actual_w(5);
results_P(7:end,cs)=errors_w(P_select,1);
results_PQ([1 3:6],cs)=actual_w([5 4 1 2 3])';
results_PQ(2,cs)=rp_exact-actual_w(5);
results_PQ(7:end,cs)=errors_w(P_select,5);
results_PQ(10,cs)=results_P(10,cs);
results_Q(1,cs)=rp_exact-actual_w(5);
results_Q(2:end,cs)=errors_w(Q_select);

results_P_export=compose('%1.3g',results_P);
results_Q_export=compose('%1.3g',results_Q);
results_PQ_export=compose('%1.3g',results_PQ);

results_short=[results_PQ(9,1) results_PQ(11,1) results_PQ(2,2:end); results_PQ(10,1) results_PQ(12,1)  results_PQ(13,2:end)];
results_short_export=compose('%1.3g',results_short);
% results=nan(5,9);
%  eig_results=nan(5,8);
%  num_results=nan(4,8);
%  impact_results=nan(9,8);
% 
% 
% 
%  results([1 3],1)=[7.8; 100*0.00566];
%  load('exact_results')
%  results([1 3],2)=[rp_exact; 100*log_c_std_exact];
%  eig_results(:,1)=[eig_seperation_exact; eigval_exact];
%  num_results(:,1)=[eig_seperation_exact; errors_exact.conditioning; errors_exact.lower_backward_error; errors_exact.upper_backward_error];
%  impact_results(:,1)=[impact_errors_exact(3,:)'; impact_errors_exact(1,:)'; impact_errors_exact(1,:)'];
%  load('solab_results')
%  results(:,3)=[rp_solab; rp_exact-rp_solab; 100*log_c_std_solab;max(max(abs(X_exact-ALPHA_ZS_solab)));max(max(abs(Q_exact-BETA_ZS_solab)))];
%  eig_results(:,2)=[eig_seperation_solab; eigval_exact-eigval_solab];
%  num_results(:,2)=[eig_seperation_solab; errors_solab.conditioning; errors_solab.lower_backward_error; errors_solab.upper_backward_error];
%   impact_results(:,2)=[impact_errors_solab(3,:)'; impact_errors_solab(1,:)'; impact_errors_solab(1,:)'];
% load('gensys_results')
%  results(:,4)=[rp_gensys; rp_exact-rp_gensys; 100*log_c_std_gensys;max(max(abs(X_exact-ALPHA_ZS_gensys)));max(max(abs(Q_exact-BETA_ZS_gensys)))];
%  eig_results(:,3)=[eig_seperation_gensys; eigval_exact-eigval_gensys];
%  num_results(:,3)=[eig_seperation_gensys; errors_gensys.conditioning; errors_gensys.lower_backward_error; errors_gensys.upper_backward_error];
%  impact_results(:,3)=[impact_errors_gensys(3,:)'; impact_errors_gensys(1,:)'; impact_errors_gensys(1,:)'];
% load('uhlig_results')
%   results(:,5)=[rp_uhlig; rp_exact-rp_uhlig; 100*log_c_std_uhlig;max(max(abs(X_exact-ALPHA_ZS_uhlig)));max(max(abs(Q_exact-BETA_ZS_uhlig)))];
% eig_results(:,4)=[eig_seperation_uhlig; eigval_exact-eigval_uhlig];
%  num_results(:,4)=[eig_seperation_uhlig; errors_uhlig.conditioning; errors_uhlig.lower_backward_error; errors_uhlig.upper_backward_error];
%  impact_results(:,4)=[impact_errors_uhlig(3,:)'; impact_errors_uhlig(1,:)'; impact_errors_uhlig(1,:)'];
% load('dynare_results')
%   results(:,6)=[rp_dynare; rp_exact-rp_dynare; 100*log_c_std_dynare;max(max(abs(X_exact-ALPHA_ZS_dynare)));max(max(abs(Q_exact-BETA_ZS_dynare)))];
% eig_results(:,5)=[eig_seperation_dynare; eigval_exact-eigval_dynare];
%  num_results(:,5)=[eig_seperation_dynare; errors_dynare.conditioning; errors_dynare.lower_backward_error; errors_dynare.upper_backward_error];
%  impact_results(:,5)=[impact_errors_dynare(3,:)'; impact_errors_dynare(1,:)'; impact_errors_dynare(1,:)'];
% load('aim_results')
%   results(:,7)=[rp_aim; rp_exact-rp_aim;100*log_c_std_aim;max(max(abs(X_exact-ALPHA_ZS_aim)));max(max(abs(Q_exact-BETA_ZS_aim)))];
% eig_results(:,6)=[eig_seperation_aim; eigval_exact-eigval_aim];
%  num_results(:,6)=[eig_seperation_aim; errors_aim.conditioning; errors_aim.lower_backward_error; errors_aim.upper_backward_error];
%  impact_results(:,6)=[impact_errors_aim(3,:)'; impact_errors_aim(1,:)'; impact_errors_aim(1,:)'];
% load('bp_results') 
% results(:,8)=[rp_bp; rp_exact-rp_bp;100*log_c_std_bp;max(max(abs(X_exact-ALPHA_ZS_bp)));max(max(abs(Q_exact-BETA_ZS_bp)))];
% eig_results(:,7)=[eig_seperation_bp; eigval_exact-eigval_bp];
%  num_results(:,7)=[eig_seperation_bp; errors_bp.conditioning; errors_bp.lower_backward_error; errors_bp.upper_backward_error];
% impact_results(:,7)=[impact_errors_bp(3,:)'; impact_errors_bp(1,:)'; impact_errors_bp(1,:)'];
% load('dynare_cr_results')
%   results(:,9)=[rp_dynare_cr; rp_exact-rp_dynare_cr; 100*log_c_std_dynare_cr;max(max(abs(X_exact-ALPHA_ZS_dynare_cr)));max(max(abs(Q_exact-BETA_ZS_dynare_cr)))];
% eig_results(:,8)=[eig_seperation_dynare_cr; eigval_exact-eigval_dynare_cr];
%  num_results(:,8)=[eig_seperation_dynare_cr; errors_dynare_cr.conditioning; errors_dynare_cr.lower_backward_error; errors_dynare_cr.upper_backward_error];
%  impact_results(:,8)=[impact_errors_dynare_cr(3,:)'; impact_errors_dynare_cr(1,:)'; impact_errors_dynare_cr(1,:)'];
% 
%  load('Y')
%  Y(1)=Y_1(1);Y(2)=0.99;Y(3)=0.025;Y(4)=0.36;Y(5)=Y_1(2);Y(6)=0.95;Y(7)=Y_1(3);
% %final_results=[results(1:3,:);NaN eig_results(1,:); NaN NaN max(abs(eig_results(2:4,2:end)));results(4:5,:)];
% final_results=[results(1:3,:);results(4:5,:);NaN eig_results(1,:); NaN NaN max(abs(eig_results(2:4,2:end))); NaN(4,1) num_results];
% final_results=[results(1:5,:);NaN NaN max(abs(eig_results(2:4,2:end))); NaN(4,1) num_results];
% short_results=nan(2,9);
% short_results(1,1:2)=num_results(1:2,1)';
% short_results(:,3:end)=[ results(2,3:end);num_results(end,2:end)];
%  save('tables','results','eig_results','Y','final_results','num_results')
%  disp(compose('%1.3g',final_results))
%  one_vector=[1 1 0 1 0 1 0];
%  disp(compose('%1.3e',Y'-one_vector))
%  disp(compose('%1.3g',short_results))
% 
% P_results=final_results([1:4 6:10],:);
% Q_results=[final_results([1:3 5],:);  nan(9,1) impact_results];
% disp('P results')
%  disp(compose('%1.3g',P_results))
%  disp('Q results')
%   disp(compose('%1.3g',Q_results))
%   
%   
%  save('tables','results','eig_results','Y','final_results','num_results','P_results','Q_results')