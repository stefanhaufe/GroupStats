% Multi Subject - Simulation 
% performs different two tailed tests for a mean difference in two classes
% and generates a plot with the true mean difference on the x-axis and the 
% H_O rejection rate on the y-axis. 

% Copyright (c) 2018 Irene Dowding and Stefan Haufe

%allocate storage space
ps_stoufer = zeros(N_rep, length(mu_diff)); 
ps_meanDiff_FE = zeros(N_rep, length(mu_diff)); 
ps_meanDiff_RE = zeros(N_rep, length(mu_diff)); 
ps_meanDiff_equal_FE = zeros(N_rep, length(mu_diff)); 
ps_meanDiff_equal_RE = zeros(N_rep, length(mu_diff)); 
ps_auc_FE = zeros(N_rep, length(mu_diff)); 
ps_auc_RE = zeros(N_rep, length(mu_diff)); 
ps_pooling = zeros(N_rep, length(mu_diff)); 
ps_pairedt = zeros(N_rep, length(mu_diff));

ps_t = zeros(1, N_vp); aucs = zeros(1, N_vp); auc_vars = zeros(1, N_vp); 
data_1 = cell(1, N_vp);

for idx_nrep = 1:N_rep
    fprintf('.')
    
%     generate_data;
      
    % or load pre-generated data
    load(['data/fig' figname '_data_rep' num2str(idx_nrep)])
    
    %% Apply different tests
    for idx_mu = 1:length(mu_diff) %for each true mean difference
        %add mean difference to create data_1, the final data set of class 1 
        for idx_vp = 1:N_vp
            data_1{idx_vp} = data_1_tmp{idx_vp} + mu_diffs(idx_vp, idx_mu); 
        end 
        
        %% Stouffer for t-test 
        %perform unpaired, one-tailed t-tests 
        for idx_vp = 1:N_vp          
            [h, ps_t(idx_vp)] = ttest2(data_0{idx_vp}, data_1{idx_vp}, 0.05, 'left');
        end

        %make sure nothing is exactly 0 or 1 
        ps_t(ps_t == 1) = 1 -eps;
        ps_t(ps_t == 0) = eps;

        %Stouffer test
        z_t = sum(norminv(ps_t)) / sqrt(N_vp); 
        ps_stoufer(idx_nrep, idx_mu) = 2*min(normcdf(z_t), 1 - normcdf(z_t)); 
        
        %% Inverse variance weighting of mean difference
        ds = cellfun(@mean, data_0) - cellfun(@mean, data_1);  % mean difference of both classes per subject
        v_0 = cellfun(@var, data_0);  v_1 = cellfun(@var, data_1);  % variances of both classes per subject
        vars = v_0 ./ N_epochs(:,1)'  + v_1 ./ N_epochs(:,2)'; 
        
        [p_FE, p_RE] = invvar_weighting(ds, vars);
        ps_meanDiff_FE(idx_nrep, idx_mu) = p_FE; %fixed effects
        ps_meanDiff_RE(idx_nrep, idx_mu) = p_RE; %random effects 
        
        %% Equal weighting of mean difference 
        [p_FE, p_RE] = equal_weighting(ds, vars);
        ps_meanDiff_equal_FE(idx_nrep, idx_mu) = p_FE;  
        ps_meanDiff_equal_RE(idx_nrep, idx_mu) = p_RE; 
                         
        %% Inverse variance weighting of AUC
        % AUC per subject and its variance
        for idx_vp = 1:N_vp          
            [foo, foo, foo, aucs(idx_vp)] = ... 
                perfcurve( [ones(N_epochs(idx_vp, 1), 1); zeros(N_epochs(idx_vp, 2),1)], ...
                                [data_0{idx_vp}, data_1{idx_vp}]', 1); 
             %variance using Hanley methods
             A = aucs(idx_vp); 
             Q1 = A/(2-A); 
             Q2 = 2*A^2/(1+A);
             n1 = N_epochs(idx_vp, 1); n2 =  N_epochs(idx_vp, 2);            
             auc_vars(idx_vp) = (A*(1-A) + (n1-1)*(Q1-A^2) +(n2-1)*(Q2 - A^2))/(n1*n2); 
        end 
        
        [p_FE, p_RE] = invvar_weighting(aucs-0.5, auc_vars);
        ps_auc_FE(idx_nrep, idx_mu) = p_FE; 
        ps_auc_RE(idx_nrep, idx_mu) = p_RE; 

        %% Pooling (putting all subjects together)
        [h, ps_pooling(idx_nrep, idx_mu)] = ttest2(cell2mat(data_0), cell2mat(data_1)); 

        %% Paired t-test on the mean per subject 
        [h, ps_pairedt(idx_nrep, idx_mu), foo, stats2] = ttest(cellfun(@mean, data_0), cellfun(@mean, data_1) ); 
    end
end

%% Generate Summary Plot
linewidth = 2; 
legend_strings = {};         

figure 
plot(mu_diff, mean(ps_stoufer < 0.05), 'k+-.', 'LineWidth', linewidth)
legend_strings = {legend_strings{:}, 'Stouffer'};

box on
hold on

plot(mu_diff, mean(ps_meanDiff_FE < 0.05), 'x-.', 'LineWidth', linewidth,'Color', [0.8500    0.3250    0.0980])
legend_strings = {legend_strings{:}, 'InvVar: Mean Difference FE'};

plot(mu_diff, mean(ps_meanDiff_RE < 0.05), 'x-', 'LineWidth', linewidth,'Color',  [0.8500    0.3250    0.0980] )
legend_strings = {legend_strings{:}, 'InvVar: Mean Difference RE'};

plot(mu_diff, mean(ps_meanDiff_equal_FE < 0.05), 'o-.', 'LineWidth', linewidth,'Color', [0.9290    0.6940    0.1250])
legend_strings = {legend_strings{:}, 'EqualVar: Mean Difference FE'};

plot(mu_diff, mean(ps_meanDiff_equal_RE < 0.05), 'o-', 'LineWidth', linewidth, 'Color', [0.9290    0.6940    0.1250])
legend_strings = {legend_strings{:}, 'EqualVar: Mean Difference RE'};

plot(mu_diff, mean(ps_auc_FE < 0.05), 'h-.', 'LineWidth', linewidth,'Color', [0.3    0.75    0.93])
legend_strings = {legend_strings{:}, 'InvVar: AUC FE'};

plot(mu_diff, mean(ps_auc_RE < 0.05), 'h-', 'LineWidth', linewidth,'Color', [0.3    0.75    0.93])
legend_strings = {legend_strings{:}, 'InvVar: AUC RE'};

plot(mu_diff, mean(ps_pooling < 0.05), 'v:', 'LineWidth', linewidth, 'Color',   [0    0.45    0.74])
legend_strings = {legend_strings{:}, 'Pooling'};

plot(mu_diff, mean(ps_pairedt < 0.05), 's-', 'LineWidth', linewidth, 'Color',   [0.4660    0.6740    0.1880])
legend_strings = {legend_strings{:}, 'Naive (paired t-test)'};

set(gca, 'FontSize', 16)
xlabel('True mean difference')
ylabel('H_0 rejection rate (\alpha = 0.05)')
legend(legend_strings, 'Location', 'Southeast')
grid on
ylim([0, 1])