% load EEG data
load data/fig4_data

% proccess data of each subject
pool1 = [];
pool2 = [];
for ivp = 1:18

  % indices of emergency braking and normal driving epochs
  in1 = find(erp_r{ivp}.y(1, :));
  in2 = find(erp_r{ivp}.y(2, :));

  % area under the ROC curve related to distinguishing the two
  % epoch types, evaluated for the two channels and each time point after
  % triggering of the emergency situation
  for ii = 1:size(erp_r{ivp}.x, 1)
    for jj = 1:size(erp_r{ivp}.x, 2)
      [~, ~, ~, me(ii, jj, ivp)] = ... 
                perfcurve( vec(erp_r{ivp}.y(1, :)), vec(erp_r{ivp}.x(ii, jj, :)), 1);
    end
  end
  
  n1 = length(in1); n2 =  length(in2);     
  
  % variance of the AUC according to Hanley 1982
  A = me(:, :, ivp); 
  Q1 = A./(2-A); 
  Q2 = 2*A.^2./(1+A);
  va(:, :, ivp) = (A.*(1-A) + (n1-1)*(Q1-A.^2) +(n2-1)*(Q2 - A.^2))./(n1*n2); 
  
  % means and variances epoch type 1 for subject-wise parametric test (AUC)
  me1(:, :, ivp) = mean(erp_r{ivp}.x(:, :, in1), 3);
  va1(:, :, ivp) = var(erp_r{ivp}.x(:, :, in1), [], 3)/size(erp_r{ivp}.x(:, :, in1), 3);
  
  % means and variances epoch type 2 for subject-wise parametric test (mean
  % difference)
  me2(:, :, ivp) = mean(erp_r{ivp}.x(:, :, in2), 3);
  va2(:, :, ivp) = var(erp_r{ivp}.x(:, :, in2), [], 3)/size(erp_r{ivp}.x(:, :, in2), 3);
  
  % collect pooled data
  pool1 = cat(3, pool1, erp_r{ivp}.x(:, :, in1));
  pool2 = cat(3, pool2, erp_r{ivp}.x(:, :, in2));
end

% fixed and random effects analysis for subject-wise non-parametric statistic (AUC)
me = me - 0.5; 
for ii = 1:size(me, 1)
  for jj = 1:size(me, 2)
    [p_FE(ii, jj), p_RE(ii, jj), z_FE(ii, jj), z_RE(ii, jj)] = invvar_weighting(me(ii, jj, :), va(ii, jj, :));
    [p_FE_eq(ii, jj), p_RE_eq(ii, jj), z_FE_eq(ii, jj), z_RE_eq(ii, jj)] = equal_weighting(me(ii, jj, :), va(ii, jj, :));
  end
end


% fixed and random effects analysis for subject-wise mean difference (parametric)  
me_t = me1-me2;
va_t = va1+va2;
for ii = 1:size(me, 1)
  for jj = 1:size(me, 2)
    [p_FE_t(ii, jj), p_RE_t(ii, jj), z_FE_t(ii, jj), z_RE_t(ii, jj)] = invvar_weighting(me_t(ii, jj, :), va_t(ii, jj, :));
  end
end

% naive t-test on subject means without considering variances
[h, p, ci, stats] = ttest(me1, me2, 'dim', 3);
z_naive = norminv(t_cdf(stats.tstat, stats.df(1)), 0, 1);

% colors
ColOrd = get(0,'DefaultAxesColorOrder');
set(0,'DefaultAxesColorOrder',ColOrd([2 1 3:7], :))

% EEG vs EMG discriminability based on parametric subject-level tests and
% inverse-variance weighting
figure(1);clf;
subplot(2, 2, 1)
plot(erp_r{1}.t, abs(z_RE_t), 'linewidth', 2);
legend('EEG', 'EMG', 'Location','NorthWest')
hold on
ax = gca;
% ax.ColorOrderIndex = 1;
% plot(erp_r{1}.t, abs(z_RE), 'linewidth', 2);
grid on
set(gca, 'fontsize', 14)
ylabel('|z|')
xlabel('ms')
ylim([0 20])
xlim([0 800])
title('EEG vs. EMG')

% random vs. fixed effect analysis using parametric subject-level tests and
% inverse-variance weighting
subplot(2, 2, 2)
plot(erp_r{1}.t, abs(z_RE_t(:, 1)), 'linewidth', 2);
hold on
ax = gca;
ax.ColorOrderIndex = 3;
plot(erp_r{1}.t, abs(z_FE_t(:, 1)), 'linewidth', 2);
legend('Random', 'Fixed', 'Location','NorthWest')
grid on
set(gca, 'fontsize', 14)
ylabel('|z|')
xlabel('ms')
ylim([0 40])
xlim([0 800])
title('Random vs. Fixed Effects')

% random effect analysis using parametric subject-level tests and
% inverse-variance weighting VS. naive t-test on subject-means
subplot(2, 2, 3)
plot(erp_r{1}.t, abs(z_RE_t(:, 1)), 'linewidth', 2);
hold on
ax = gca;
ax.ColorOrderIndex = 5;
plot(erp_r{1}.t, abs(z_naive(:, 1)), 'linewidth', 2);
legend('Inverse Variance', 'Naive', 'Location','NorthWest')
grid on
set(gca, 'fontsize', 14)
ylabel('|z|')
xlabel('ms')
ylim([0 10])
xlim([0 800])
title('InvVar vs. Naive Summary Stat')

% parametric (mean-diff) VS. nonparametric (AUC) subject-level test using
% random effect analysis using inverse-variance weighting 
subplot(2, 2, 4)
plot(erp_r{1}.t, abs(z_RE_t(:, 1)), 'linewidth', 2);
hold on
ax = gca;
ax.ColorOrderIndex = 6;
plot(erp_r{1}.t, abs(z_RE(:, 1)), 'linewidth', 2);
legend('Mean Diff', 'AUC', 'Location','NorthWest')
grid on
set(gca, 'fontsize', 14)
ylabel('|z|')
xlabel('ms')
ylim([0 10])
xlim([0 800])
title('Gaussian vs. Nonparametric')


% export_fig(['erp1'], '-r300', '-a2'); 



