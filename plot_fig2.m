% Copyright (c) 2018 Irene Dowding and Stefan Haufe

clear all; close all;

% fig 2A (Fixed Effects, S=5)
N_rep = 100; %number of repetitions of the simulation
N_vp = 5; %number of subjects
N_epochs_range = [50 80]; % range of number of sample size per subject and class (uniformly drawn)
mu_diff = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3];  %mean difference of two classes 
                                                 %(x-axis of the generated plot)
normal = 1; %if 1, data is from a normal distribution, otherwise from an F-distribution
sigma_epochs_range = [0.5 2]; %each subject s may have a different st.d. sigma_s around their class means. 
                     %The parameters sets the range of sigma_s per subject 
                     % e.g. [1 1] - each subject has sigma_s 1
                     % e.g. [0.5 2] each subjects sigma_s is uniformy drawn from [0.5, 2] 
vp_mu_range = 6; %each subject may have a different mean (average for both classes)
                 %if set to 0, each subject has the same mean
sigma_rand = 0; %st.d. of mu_diff per subject 
                  %if 0 - fixed effect model, each subject has the same mean difference.
var_corr = 0; % correlation of variance with effect size (mean difference)
              %  0 - null              
              %  1 - positive
              % -1 - negative    
              
figname = '2A';
main_simulation
title('Fixed Effect Model, S = 5')
export_fig(['figures/fig2A'], '-r300', '-a2'); 

% fig 2B (Fixed Effects, S=20)
N_vp = 20; %number of subjects
figname = '2B';
main_simulation
title('Fixed Effect Model, S = 20')
export_fig(['figures/fig2B'], '-r300', '-a2'); 

% fig 2D (Random Effects, S=20)
sigma_rand = 0.2;
figname = '2D';
main_simulation
title('Random Effects Model, S = 20')
export_fig(['figures/fig2C'], '-r300', '-a2'); 

% fig 2C (Random Effects, S=5)
N_vp = 20; %number of subjects
figname = '2C';
main_simulation
title('Random Effects Model, S = 5')
export_fig(['figures/fig2D'], '-r300', '-a2'); 
