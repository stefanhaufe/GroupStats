function [p_FE, p_RE, z_FE, z_RE] = invvar_weighting(effects, vars)
% Compute Fixed Effects (FE) and Random Effects (RE) p-values using 
% inverse variance weighting
% IN:  effects      - array of effect sizes (one per subject)
%      vars         - array of estimated variance of effect size (one per subject)
% OUT: p_FE, p_RE   - p value of H_0 : mean effect == 0 , for fixed and
%                       random effects analysis
%
% Copyright (c) 2018 Irene Dowding and Stefan Haufe

effects = effects(:); 
vars = vars(:); 

%% Compute p-value of fixed effects analysis
weights = 1./vars; 
meaneffect = sum(weights.* effects)/sum(weights);  %weighted mean effect
meanvariance = 1./sum(weights); %variance of the weighted mean
z_FE = meaneffect / sqrt(meanvariance); 
p_FE = 2*min(normcdf(z_FE), 1 - normcdf(z_FE)); 

%% Compute p-value of random effects analysis
sigma2_rand = getBetweenSubjectVarianceWithDerSimonianAndLaird(effects, vars); 
weights = 1./(vars + sigma2_rand); 
meaneffect = sum(weights.* effects)/sum(weights); %weighted mean effect
meanvariance = 1./sum(weights); %variance of the weighted mean
z_RE = meaneffect / sqrt(meanvariance); 
p_RE = 2*min(normcdf(z_RE), 1 - normcdf(z_RE)); 

function sigma2_rand = getBetweenSubjectVarianceWithDerSimonianAndLaird(effects, vars)
    % Compute random effect variance using DerSimonian and Laird
    N_subjects = length(effects); 
    a = 1./vars;
    Q = a' * (effects - sum(a.* effects)/sum(a)).^2;  
    sigma2_rand = max(0, (Q-N_subjects+1)/(sum(a) - sum(a.^2)/sum(a))); 

    


