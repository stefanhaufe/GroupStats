function [p_FE, p_RE, z_FE, z_RE] = equal_weighting(effects, vars)
% Compute Fixed Effects (FE) and Random Effects (RE) p-values using 
% equal weighting of the effect sizes
% IN:  effects      - array of effect sizes (one per subject)
%      vars         - array of estimated variance of effect size (one per subject)
% OUT: p_FE, p_RE   - p value of H_0 : mean effect == 0 , for fixed and
%                       random effects analysis
%
% Copyright (c) 2018 Irene Dowding and Stefan Haufe


effects = effects(:); 
vars = vars(:); 
N_vp = length(effects);

%% Compute p-value according to fixed effects analysis
weights = repmat(1/N_vp, N_vp, 1);
meaneffect = sum(weights.* effects)/sum(weights);  %weighted mean effect
z_FE = meaneffect / sqrt(sum(weights.^2.*vars)); 
p_FE = 2*min(normcdf(z_FE), 1 - normcdf(z_FE)); 

%% Compute p-value of random effects analysis
sigma2_rand = getBetweenSubjectVarianceWithDerSimonianAndLaird(effects, vars); 
z_RE = meaneffect / sqrt(sum(weights.^2.*(vars + sigma2_rand))); 
p_RE = 2*min(normcdf(z_RE), 1 - normcdf(z_RE)); 
    
function sigma2_rand = getBetweenSubjectVarianceWithDerSimonianAndLaird(effects, vars)
    % returns random effect variance computed using DerSimonian and Laird
    N_subjects = length(effects); 
    a = 1./vars;
    Q = a' * (effects - sum(a.* effects)/sum(a)).^2;  
    sigma2_rand = max(0, (Q-N_subjects+1)/(sum(a) - sum(a.^2)/sum(a))); 