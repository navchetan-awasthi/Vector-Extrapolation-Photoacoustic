%% This generates the results for the experimental phantom

A_b=A_bfTriangle;
b=pt_data(:);
k=2;
L=100;
tol=10;


x0=zeros(40000,1);
i=4;
[x1_mpe_bv{i}, residuals, final_bv(i)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'MPE',tol);
 figure;imshow(reshape(x1_mpe_bv{i},200,200),[]);
[x1_rre_bv{i+1}, residuals, final_bv(i+1)] = my_extrapolate_rsd(A_b, b, x0, k, L, 'RRE',tol);
 figure;imshow(reshape(x1_rre_bv{i+1},200,200),[]);
 x0=A_b'*b;
[x1_sd_bv{i+2},residuals,final_bv(i+2)]=rsd(A_b,b,x0,tol);
figure;imshow(reshape(x1_sd_bv{i+2},200,200),[]);


x0=zeros(40000,1);
i=7;
[x1_mpe_bv{i}, residuals, final_bv(i)] = my_extrapolate_tv_hair(V,S,A_b, b, x0, k, L, 'MPE',tol);
 figure;imshow(reshape(x1_mpe_bv{i},200,200),[])
[x1_rre_bv{i+1}, residuals, final_bv(i+1)] = my_extrapolate_tv_hair(V,S,A_b, b, x0, k, L, 'RRE',tol);
 figure;imshow(reshape(x1_rre_bv{i+1},200,200),[])
[x1_sd_bv{i+2},final_bv(i+2)]=tv_hair(V,S,b,A_b,tol)
figure;imshow(reshape(x1_sd_bv{i+2},200,200),[])



i=7;
x_img = x1_mpe_bv{i};
imgmax = max(x_img(:));
imgmin = min(x_img(:));
nse = std(x_img(:));
snr = 20*log((imgmax-imgmin)/nse);
x_img = x1_rre_bv{i+1};
imgmax = max(x_img(:));
imgmin = min(x_img(:));
nse = std(x_img(:));
snr = 20*log((imgmax-imgmin)/nse);
x_img = x1_sd_bv{i+2};
imgmax = max(x_img(:));
imgmin = min(x_img(:));
nse = std(x_img(:));
snr = 20*log((imgmax-imgmin)/nse);


