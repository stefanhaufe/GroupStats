
mu_diffs = zeros(N_vp, length(mu_diff)); %mean difference per subject 
data_0 = cell(1, N_vp); data_1_tmp = cell(1, N_vp);  

%% Generate Data 
%determine random number of epochs per subject and class
% Copyright (c) 2018 Irene Winkler and Stefan Haufe

N_epochs = randi(N_epochs_range, [N_vp, 2]);

%generate data for class 0 and class 1 
for idx_vp = 1:N_vp
    if normal == 1
        data_0{idx_vp} = randn(1, N_epochs(idx_vp, 1)); 
        data_1_tmp{idx_vp} = randn(1, N_epochs(idx_vp, 2));
    else
        %F-distribution
        d1 = 2; d2 = 5; 
        data_0{idx_vp}  = frnd(d1, d2, [1, N_epochs(idx_vp, 1)]) - d2/(d2-2); 
        data_1_tmp{idx_vp}  = frnd(d1, d2, [1, N_epochs(idx_vp, 2)]) - d2/(d2-2);
        s = sqrt(2*d2^2*(d1+d2-2)/(d1*(d2-2)^2*(d2-4))); 
        data_0{idx_vp}  = 1/s*data_0{idx_vp}; %make it std. 1
        data_1_tmp{idx_vp}  = 1/s*data_1_tmp{idx_vp} ; %make it std. 1
    end

    %standard deviation for each subject
    s = sigma_epochs_range(1) + rand(1)* (sigma_epochs_range(2) - sigma_epochs_range(1));
    data_0{idx_vp} = s * data_0{idx_vp};
    data_1_tmp{idx_vp} = s * data_1_tmp{idx_vp}; 

    %mean for each subject 
    mu = rand(1)*vp_mu_range - vp_mu_range/2; 
    data_0{idx_vp} = data_0{idx_vp} + mu;
    data_1_tmp{idx_vp} = data_1_tmp{idx_vp} + mu;

    %mean differences for each subject (same for fixed effect model)
    mu_diffs(idx_vp, :) = mu_diff + randn(1, length(mu_diff))* sigma_rand; 
end
