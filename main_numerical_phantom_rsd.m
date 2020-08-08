%% This genertes the results for Numerical Phantom for Regularized Steepest Descent Methods

%% 20 db simulations

load detectors_20_2(100).mat
k=2;
L=100;
tol=10;
i=1;
x0=A_b'*sdn2_v_bv;
b=sdn2_v_bv;
BV2=double(BV2_bv);
[x1_mpe_bv{i}, residuals_bv{i}, final_bv(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_bv(i),cnr_bv(i)]=comp(BV2,x1_mpe_bv{i});
[RMSE_bv(i)] = my_rmse(reshape(x1_mpe_bv{i},201,201),BV2);
[x1_rre_bv{i+1}, residuals_bv{i+1}, final_bv(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_bv(i+1),cnr_bv(i+1)]=comp(BV2,x1_rre_bv{i+1});
[RMSE_bv(i+1)] = my_rmse(reshape(x1_rre_bv{i+1},201,201),BV2);
[x1_sd_bv{i+2},residuals_bv{i+2},final_bv(i+2)]=rsd(A_b,b,x0,tol);
[pc_bv(i+2),cnr_bv(i+2)]=comp(BV2,x1_sd_bv{i+2});
[RMSE_bv(i+2)] = my_rmse(reshape(x1_sd_bv{i+2},201,201),BV2);

x0=A_b'*sdn2_v_der;
b=sdn2_v_der;
BV2=double(BV2_der);
[x1_mpe_der{i}, residuals_der{i}, final_der(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_der(i),cnr_der(i)]=comp(BV2,x1_mpe_der{i});
[RMSE_der(i)] = my_rmse(reshape(x1_mpe_der{i},201,201),BV2);
[x1_rre_der{i+1}, residuals_der{i+1}, final_der(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_der(i+1),cnr_der(i+1)]=comp(BV2,x1_rre_der{i+1});
[RMSE_der(i+1)] = my_rmse(reshape(x1_rre_der{i+1},201,201),BV2);
[x1_sd_der{i+2},residuals_der{i+2},final_der(i+2)]=rsd(A_b,b,x0,tol);
[pc_der(i+2),cnr_der(i+2)]=comp(BV2,x1_sd_der{i+2});
[RMSE_der(i+2)] = my_rmse(reshape(x1_sd_der{i+2},201,201),BV2);


x0=A_b'*sdn2_v_pat;
b=sdn2_v_pat;
BV2=double(BV2_pat);
[x1_mpe_pat{i}, residuals_pat{i}, final_pat(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_pat(i),cnr_pat(i)]=comp(BV2,x1_mpe_pat{i});
[RMSE_pat(i)] = my_rmse(reshape(x1_mpe_pat{i},201,201),BV2);
[x1_rre_pat{i+1}, residuals_pat{i+1}, final_pat(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_pat(i+1),cnr_pat(i+1)]=comp(BV2,x1_rre_pat{i+1});
[RMSE_pat(i+1)] = my_rmse(reshape(x1_rre_pat{i+1},201,201),BV2);
[x1_sd_pat{i+2},residuals_pat{i+2},final_pat(i+2)]=rsd(A_b,b,x0,tol);
[pc_pat(i+2),cnr_pat(i+2)]=comp(BV2,x1_sd_pat{i+2});
[RMSE_pat(i+2)] = my_rmse(reshape(x1_sd_pat{i+2},201,201),BV2);



%% 40 db simulations
load detectors_40_2(100).mat
k=2;
L=100;
tol=1;
i=4;
x0=A_b'*sdn2_v_bv;
b=sdn2_v_bv;
BV2=double(BV2_bv);
[x1_mpe_bv{i}, residuals_bv{i}, final_bv(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_bv(i),cnr_bv(i)]=comp(BV2,x1_mpe_bv{i});
[RMSE_bv(i)] = my_rmse(reshape(x1_mpe_bv{i},201,201),BV2);
[x1_rre_bv{i+1}, residuals_bv{i+1}, final_bv(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_bv(i+1),cnr_bv(i+1)]=comp(BV2,x1_rre_bv{i+1});
[RMSE_bv(i+1)] = my_rmse(reshape(x1_rre_bv{i+1},201,201),BV2);
[x1_sd_bv{i+2},residuals_bv{i+2},final_bv(i+2)]=rsd(A_b,b,x0,tol);
[pc_bv(i+2),cnr_bv(i+2)]=comp(BV2,x1_sd_bv{i+2});
[RMSE_bv(i+2)] = my_rmse(reshape(x1_sd_bv{i+2},201,201),BV2);

x0=A_b'*sdn2_v_der;
b=sdn2_v_der;
BV2=double(BV2_der);
[x1_mpe_der{i}, residuals_der{i}, final_der(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_der(i),cnr_der(i)]=comp(BV2,x1_mpe_der{i});
[RMSE_der(i)] = my_rmse(reshape(x1_mpe_der{i},201,201),BV2);
[x1_rre_der{i+1}, residuals_der{i+1}, final_der(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_der(i+1),cnr_der(i+1)]=comp(BV2,x1_rre_der{i+1});
[RMSE_der(i+1)] = my_rmse(reshape(x1_rre_der{i+1},201,201),BV2);
[x1_sd_der{i+2},residuals_der{i+2},final_der(i+2)]=rsd(A_b,b,x0,tol);
[pc_der(i+2),cnr_der(i+2)]=comp(BV2,x1_sd_der{i+2});
[RMSE_der(i+2)] = my_rmse(reshape(x1_sd_der{i+2},201,201),BV2);


x0=A_b'*sdn2_v_pat;
b=sdn2_v_pat;
BV2=double(BV2_pat);
[x1_mpe_pat{i}, residuals_pat{i}, final_pat(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_pat(i),cnr_pat(i)]=comp(BV2,x1_mpe_pat{i});
[RMSE_pat(i)] = my_rmse(reshape(x1_mpe_pat{i},201,201),BV2);
[x1_rre_pat{i+1}, residuals_pat{i+1}, final_pat(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_pat(i+1),cnr_pat(i+1)]=comp(BV2,x1_rre_pat{i+1});
[RMSE_pat(i+1)] = my_rmse(reshape(x1_rre_pat{i+1},201,201),BV2);
[x1_sd_pat{i+2},residuals_pat{i+2},final_pat(i+2)]=rsd(A_b,b,x0,tol);
[pc_pat(i+2),cnr_pat(i+2)]=comp(BV2,x1_sd_pat{i+2});
[RMSE_pat(i+2)] = my_rmse(reshape(x1_sd_pat{i+2},201,201),BV2);


%% 60 db simulations

load detectors_60_2(100).mat
k=2;
L=100;
tol=0.1;
i=7;
x0=A_b'*sdn2_v_bv;
b=sdn2_v_bv;
BV2=double(BV2_bv);
[x1_mpe_bv{i}, residuals_bv{i}, final_bv(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_bv(i),cnr_bv(i)]=comp(BV2,x1_mpe_bv{i});
[RMSE_bv(i)] = my_rmse(reshape(x1_mpe_bv{i},201,201),BV2);
[x1_rre_bv{i+1}, residuals_bv{i+1}, final_bv(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_bv(i+1),cnr_bv(i+1)]=comp(BV2,x1_rre_bv{i+1});
[RMSE_bv(i+1)] = my_rmse(reshape(x1_rre_bv{i+1},201,201),BV2);
[x1_sd_bv{i+2},residuals_bv{i+2},final_bv(i+2)]=rsd(A_b,b,x0,tol);
[pc_bv(i+2),cnr_bv(i+2)]=comp(BV2,x1_sd_bv{i+2});
[RMSE_bv(i+2)] = my_rmse(reshape(x1_sd_bv{i+2},201,201),BV2);

x0=A_b'*sdn2_v_der;
b=sdn2_v_der;
BV2=double(BV2_der);
[x1_mpe_der{i}, residuals_der{i}, final_der(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_der(i),cnr_der(i)]=comp(BV2,x1_mpe_der{i});
[RMSE_der(i)] = my_rmse(reshape(x1_mpe_der{i},201,201),BV2);
[x1_rre_der{i+1}, residuals_der{i+1}, final_der(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_der(i+1),cnr_der(i+1)]=comp(BV2,x1_rre_der{i+1});
[RMSE_der(i+1)] = my_rmse(reshape(x1_rre_der{i+1},201,201),BV2);
[x1_sd_der{i+2},residuals_der{i+2},final_der(i+2)]=rsd(A_b,b,x0,tol);
[pc_der(i+2),cnr_der(i+2)]=comp(BV2,x1_sd_der{i+2});
[RMSE_der(i+2)] = my_rmse(reshape(x1_sd_der{i+2},201,201),BV2);


x0=A_b'*sdn2_v_pat;
b=sdn2_v_pat;
BV2=double(BV2_pat);
[x1_mpe_pat{i}, residuals_pat{i}, final_pat(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
[pc_pat(i),cnr_pat(i)]=comp(BV2,x1_mpe_pat{i});
[RMSE_pat(i)] = my_rmse(reshape(x1_mpe_pat{i},201,201),BV2);
[x1_rre_pat{i+1}, residuals_pat{i+1}, final_pat(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
[pc_pat(i+1),cnr_pat(i+1)]=comp(BV2,x1_rre_pat{i+1});
[RMSE_pat(i+1)] = my_rmse(reshape(x1_rre_pat{i+1},201,201),BV2);
[x1_sd_pat{i+2},residuals_pat{i+2},final_pat(i+2)]=rsd(A_b,b,x0,tol);
[pc_pat(i+2),cnr_pat(i+2)]=comp(BV2,x1_sd_pat{i+2});
[RMSE_pat(i+2)] = my_rmse(reshape(x1_sd_pat{i+2},201,201),BV2);

