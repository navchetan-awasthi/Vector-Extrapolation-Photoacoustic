%% Used to calculate PC and CNR
function [PC1, CNR1]=comp(BV2,x1)
x1=x1;
    x_true=double(BV2);
    BV2=x_true;
    COV =  cov(BV2,x1);
    Nr =  COV(2,1) ;
    Dr =sqrt(COV(1,1) *COV(2,2))  ;
    PC1 = Nr/Dr;
    CNR1 = CONTRAST_NOISE_RATIO(BV2,x1);
    
