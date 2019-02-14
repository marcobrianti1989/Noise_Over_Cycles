%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 - Forecasted growth rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real GDP
Delta_RGDP_t        = RGDP5_SPF./RGDP1_SPF - ones(length(RGDP1_SPF),1);
Delta_RDGP_t1       = RGDP6_SPF./RGDP2_SPF - ones(length(RGDP1_SPF),1);
% Nominal GDP
Delta_NGDP_t        = NGDP5_SPF./NGDP1_SPF - ones(length(NGDP1_SPF),1);
Delta_NDGP_t1       = NGDP6_SPF./NGDP2_SPF - ones(length(NGDP1_SPF),1);
% Real Cons
Delta_RCONS_t       = RCONS5_SPF./RCONS1_SPF - ones(length(RCONS1_SPF),1);
Delta_RCONS_t1      = RCONS6_SPF./RCONS2_SPF - ones(length(RCONS1_SPF),1);
%Industrial Production
Delta_INDPROD_t     = INDPROD5_SPF./INDPROD1_SPF - ones(length(INDPROD1_SPF),1);
Delta_INDPROD_t1    = INDPROD6_SPF./INDPROD2_SPF - ones(length(INDPROD1_SPF),1);
%Investment is the sum between residential and non residential investment
Delta_RINV_t        = (RRESINV5_SPF + RNRESIN5_SPF)./(RRESINV1_SPF + RNRESIN1_SPF)  - ones(length(RRESINV1_SPF),1);
Delta_RINV_t1       = (RRESINV6_SPF + RNRESIN6_SPF)./(RRESINV2_SPF + RNRESIN2_SPF)  - ones(length(RRESINV1_SPF),1);
% CPI
Delta_CPI_t         = CPI5_SPF;% - CPI1_SPF;
Delta_CPI_t1        = CPI6_SPF;% - CPI2_SPF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2 - Revision in forecast growth rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z1                  = [NaN; Delta_RGDP_t(2:end) - Delta_RDGP_t1(1:end-1)];
Z2                  = [NaN; Delta_NGDP_t(2:end) - Delta_NDGP_t1(1:end-1)];
Z3                  = [NaN; Delta_RCONS_t(2:end) - Delta_RCONS_t1(1:end-1)];
Z4                  = [NaN; Delta_INDPROD_t(2:end) - Delta_INDPROD_t1(1:end-1)];
Z5                  = [NaN; Delta_RINV_t(2:end) - Delta_RINV_t1(1:end-1)];
Z6                  = [NaN; Delta_CPI_t(2:end)];% - Delta_CPI_t1(1:end-1)];
Z7                  = [NaN; diff(MichIndexConfidence)];
