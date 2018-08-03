function [IRF] = stlp(Y,X)
    B               = X'*X\(X'*Y);
    IRF(h)          = B(1); %here we are fixing the probability over the IRFs - need to relax
    IRF_R(h)        = B(1) + B(2); %this is correct, it can be proved mathematically in 2 steps
    IRF_L(h)        = BL(1); %unconditional OLS, not controlling for smooth transition F(z)
    res_uncond{h}   = Y - XL*BL;
    Rsquared(h)     = 1 - var(res_uncond{h})/var(Y);
    BL_store(:,h)   = BL;
    XL_store{h}     = XL;
    SE_store(h)     = se(1);%this is only for the case of expansion -- to correct 
    SEL_store(h)    = sel(1);
end




