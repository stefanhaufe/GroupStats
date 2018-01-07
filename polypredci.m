function [ p, yhat, ci ] = polypredci( x, y, n, alfa, xv )
% POLYPREDCI Calculates the prediction intervals for a polynomial
%            regression produced by ?polyfit?.
%
% DESCRIPTION:
%   The function takes the vector data in ?x? and ?y? arguments,
%       polynomial degree ?n?, and significance level ?alfa?, then uses
%       ?polyfit? and ?polyval? to fit the data and estimte the parameters.
%       It then calculates the residuals, calculates the ?t? value using
%       the same routines as in ?tstat3? (see the ?tstat3? function for
%       details).  They are included in this function so ?tstat3? is not
%       required.  It then calculates the prediction intervals, and returns
%       the estimated parameters, fitted line, and prediction intervals.
%   NOTE: The ?prediction interval? are the confidence intervals on the
%       fitted regression line, similar to the confidence intervals on any
%       other estimated statistical parameter.  They are not confidence
%       intervals on the data used in the regression.
%
% USE:
%   INPUT ARGUMENTS:
%       x       Data independent variable vector
%       y       Data dependent variable vector
%       n       Polynomial degree
%       alfa    Significance level (greater than 0 and less than 1)
%               (Default = 0.95)
%       xv      Optional longer vector to use to calculate confidence
%               intervals with increased resolution
%
%   OUTPUTS:
%       p       Regression parameters (polynomial coefficients, ?polyfit?)
%       yhat    Fitted regression line (?polyval?)
%       ci      Prediction confidence limits (vector)
% 
%   NOTE: If ?xv? is provided, ?yhat? and ?ci? will be the same length 
%         as ?xv?.  
%         If ?xv? is not provided, ?yhat? and ?ci? will be the same length 
%         as ?x?.  
%
% SOURCE: https://en.wikipedia.org/wiki/Simple_linear_regression
% AUTHOR: Star Strider ? 2016 06 11
% REVISION HISTORY: 2016 06 11
%

% %% ?> TO DO: Trap to be sure the order of the polynomial is >=1 and & # data > order+2

x = x(:);
y = y(:);

L = length(x);                                                      % Vector Length

if L ~= length(y)
    error('\n\tVectors must be the same length.  Length x = %.0f,  Length y = %.0f\n', L, length(y))
end

if n < 1
    error('\n\tPolynomial order must be >= 1\n')
end

if L < (n+2)
    error('\n\tPolynomial order must be <= %.0f for vector length %.0f\n', L-2, L)
end

if nargin < 4
    alfa = 0.95;
end

% Check for out-of-range values for alpha and substitute if necessary:
if alfa < 1.0E-010
    alfa = 1.0E-010;
elseif alfa > (1 - 1.0E-010)
    alfa = 1 - 1.0E-010;
end

if alfa > 0.5
    prob = (1 - alfa)/(2);
elseif alfa < 0.5
    prob = alfa/(2);
end
% prob

p    = polyfit(x, y, n);                                            % Fit & Evaluate Polynomial
yhat = polyval(p, x);


% This calculates the inverse t-distribution (parameters given the
%   probability ?alpha? and degrees of freedom ?v?:
tdist2T = @(t,v) (1-betainc(v/(v+t^2), v/2, 0.5));                  % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;                              % 1-tailed t-distribution
t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist2T(tval,v)), 0);  % T-Statistic Given Probability ?alpha? & Degrees-Of-Freedom ?v?

% Calculate Confidence Intervals ?ci?
d = y(:) - yhat;                                                    % Calculate Residuals
mx = mean(x);                                                       % Mean
% L = length(x);                                                      % Vector Length
v = L-(n+1);                                                        % Degrees-Of-Freedom
tv = abs(t_inv(prob , v));                                           % t-Value
ci = tv*sqrt((sum(d.^2)/v) .* ((1/L) + ((x(:)-mx).^2)./(sum((x(:)-mx).^2))));       % Calculate Confidence Intervals With ?x?

if nargin == 5                                                      % Check For Presence Of ?xv?
    yhat = polyval(p, xv);                                          % Calculate ?yhat? With ?xv?
    ci = tv*sqrt((sum(d.^2)/v) .* ((1/L) + ((xv(:)-mx).^2)./(sum((xv(:)-mx).^2)))); % Calculate Confidence Intervals With ?xv?
end

end
% =========================== END: polypredci.m ===========================

