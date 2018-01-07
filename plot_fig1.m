% Copyright (c) 2018 Irene Winkler and Stefan Haufe

addpath external

% generate new data

% number of subjects/groups
ng = 4; 

% number of samples per group
m = 20; 

% group means
gme = 15*randn(ng, 1); 

% mean differences between two classes for each group
gdiff = reshape([-1 1], 1, 1, 2); 

% data (m samples, ng groups, two classes) 
x = 2*randn(ng, m, 2) + repmat(gme, 1, m, 2) + 1*repmat(gdiff, ng, m); 

% % or load the original data used in the paper
% load('fig1_data', 'ng', 'm', 'gme', 'gdiff', 'x');

% t-test on individual subjects vs. t-test on pooled data of all subjects
figure; 
% individual
for ii = 1:ng
  subplot(3, ng+1, ii)
  hold on
  x1 = squeeze(x(ii, :, 1));
  errorbar(1, mean(x1), std(x1)/sqrt(m), '.', 'linewidth', 5, 'markersize', 3)
  x1 = squeeze(x(ii, :, 2));
  errorbar(2, mean(x1), std(x1)/sqrt(m), '.', 'linewidth', 5, 'markersize', 3)
  xlim([0.5 2.5])
  ylim([min(x(:)) max(x(:))])
  grid on
  set(gca, 'xtick', 1:2, 'xticklabel', {'',''})
  title(['S' num2str(ii)])
end

% pooled
x2 = reshape(x, ng*m, 2);
subplot(3, ng+1, ng+1)
hold on
x2 = reshape(x(:, :, 1), [], 1);
errorbar(1, mean(x2), std(x2)/sqrt(ng*m), '.', 'linewidth', 5, 'markersize', 3)
x2 = reshape(x(:, :, 2), [], 1);
errorbar(2, mean(x2), std(x2)/sqrt(ng*m), '.', 'linewidth', 5, 'markersize', 3)
xlim([0.5 2.5])
ylim([min(x(:)) max(x(:))])
grid on
set(gca, 'xtick', 1:2, 'xticklabel', {'', ''})
title('Pooled')

% different zoom level individual
for ii = 1:ng
  subplot(3, ng+1, ng+1+ii)
  hold on
  x1 = squeeze(x(ii, :, 1));
  errorbar(1, mean(x1), std(x1)/sqrt(m), '.', 'linewidth', 5, 'markersize', 3)
  x1 = squeeze(x(ii, :, 2));
  errorbar(2, mean(x1), std(x1)/sqrt(m), '.', 'linewidth', 5, 'markersize', 3)
  xlim([0.5 2.5])
  grid on
  set(gca, 'xtick', 1:2, 'xticklabel', {['X_' num2str(ii)], ['Y_' num2str(ii)]})
  [h p] = ttest2(squeeze(x(ii, :, 1)), squeeze(x(ii, :, 2)), 'vartype', 'unequal');
  title(['p = ' num2str(p, 2)])
end

% pooled
x2 = reshape(x, ng*m, 2);
subplot(3, ng+1, 2*(ng+1))
hold on
x2 = reshape(x(:, :, 1), [], 1);
errorbar(1, mean(x2), std(x2)/sqrt(ng*m), '.', 'linewidth', 5, 'markersize', 3)
x2 = reshape(x(:, :, 2), [], 1);
errorbar(2, mean(x2), std(x2)/sqrt(ng*m), '.', 'linewidth', 5, 'markersize', 3)
xlim([0.5 2.5])
grid on
set(gca, 'xtick', 1:2, 'xticklabel', {'X', 'Y'})
[h p, ci, stats] = ttest2(reshape(x(:, :, 1), [], 1), reshape(x(:, :, 2), [], 1), 'vartype', 'unequal');
title(['p = ' num2str(p, 2)])

% export_fig(['pooling1'], '-r300', '-a2'); 

% linear regression and correlation for individual and pooled data
figure;

% plot entire data
subplot(3, ng+1, [1:ng])
plot(reshape(permute(x, [2 1 3]), [], 2), '*')
grid on
legend('X', 'Y')

% individual linear regression/correlation
for ii = 1:ng
  subplot(3, ng+1, ng+1+ii)
  hold on
  xx = squeeze(x(ii, :, 1))';
  yy = squeeze(x(ii, :, 2))';
  plot(xx, yy, '+')
  xlabel(['X_' num2str(ii)])
  if ii == 1
    ylabel(['Y_*'])
  end
  grid on
  [r p] = corr(xx, yy);
  title({['r = ' num2str(r, 2)], ['p = ' num2str(p, 2)]}')

  ax = get(gca, 'xlim');
  xval = linspace(ax(1), ax(2), 1000)';
  [p, yhat, ci ] = polypredci(xx, yy, 1, 0.95, xval); 

  plot(xval,yhat+ci,'r-.'); 
  plot(xval,yhat-ci,'r-.');
  plot(xval,yhat,'b','linewidth',2);
   
end

% same on pooled data
subplot(3, ng+1, 2*(ng+1))
hold on
xx = reshape(x(:, :, 1), [], 1);
yy = reshape(x(:, :, 2), [], 1);
plot(xx, yy, '+')
xlabel('X')
mi = min((x(:)))*1.1;
ma = max((x(:)))*1.1;
xlim([mi ma])
ylim([mi ma])
grid on
[r p] = corr(xx, yy);
title({['r = ' num2str(r, 2)], ['p = ' num2str(p, 2)]}')

ax = get(gca, 'xlim');
xval = linspace(ax(1), ax(2), 1000)';
[p, yhat, ci ] = polypredci(xx, yy, 1, 0.95, xval); 

plot(xval,yhat+ci,'r-.'); 
plot(xval,yhat-ci,'r-.');
plot(xval,yhat,'b','linewidth',2);
 
% export_fig(['pooling2'], '-r300', '-a2'); 

  