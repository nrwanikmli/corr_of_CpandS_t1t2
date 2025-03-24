% scatter_21b2015.m
% In this script, scatter plots are generated for the pair
%  - S and T_1 {-0.1430}
%  - S and T_2 {-0.0665}
%  - Cp and T_1 {-0.0245}
%  - Cp and T_2 {-0.0082}
% ... with their respective regression lines.
close all; % Before drawing, close any figures already opened
clear;     % Clear all variables
% CONFIG: Line width and font size for each curve in drawn figures
lineWidth = 2;
fontSize = 20;
% Save plots to images? Set to true if yes
saveToFile = false;
% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));
% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([% Entropy
   269.722 334.155 389.475 395.882 444.724 447.437 457.958 455.839 450.418 399.491
   499.831 513.857 508.537 507.395 506.076 512.523 500.734 520.307 509.210 513.879
   511.770 509.611 462.545 463.738 468.712 555.409 472.295 554.784 468.796 551.708
   ]', 30, 1), reshape([ % Heat capacity
   83.019 133.325 184.194 183.654 235.165 233.497 234.568 234.638 233.558 200.815
   286.182 285.056 284.037 284.088 285.148 284.595 284.870 284.503 284.785 284.740
   284.233 284.552 251.175 250.568 251.973 336.098 267.543 337.204 285.041 368.518
   ]', 30, 1),
   "S^o", "Cp" % Their labels
 };
% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
 d_f = [
   6 6 6 7 6  8 7  8 9 6 6  7  8  8  7  10 9  8  9  8  9  8  8 8 7  6  7  6  6  4
   0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8  10 12 10 12 10 10 8 8 10 12 10 20 12 16
   0 1 2 3 3  5 4  5 6 5 4  5  6  6  5  8  7  6  7  6  7  6  8 8 7  12 10 5  12 19
 ]'; % Used for computing indices based on edge endpoint degree partitions
% Define v_exp here before using it in getIndexFns 30 lower BHs (v_exp: number of vertices)
v_exp = [
    6 10 14 14 18 18 18 18 18 18 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 24 22 26 24 32
]';
% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(alpha, v_exp) | n=1:T_1^a, n=2:T_2^a, a = alpha
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (v_exp - 2)), ((5 .* v_exp - 12) ./ ((v_exp - 2) .* (v_exp - 3))), (6 ./ (v_exp - 3))].^alpha, 2); % First General Temperature Indices (T_1^a)
    @(alpha) sum(d_f .* [(4 ./ ((v_exp - 2).^2)), (6 ./ ((v_exp - 2) .* (v_exp - 3))), (9 ./ ((v_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_2^a)
}';
% Cell containing their respective labels
indexName = {"T_1" "T_2"};
% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two S and Cp
numIndices = size(getIndexFns,2); % two (T_1 & T_2)
numCurves = numData*numIndices;   % two (4 curves)
for edn = 1:numData % edn = experimental data number | 1=S, 2=Cp
  for fnn = 1:numIndices % fnn = function number | n=1:T_1, n=2:T_2
    ccFn = @(alpha) corrcoef( % Gets corrcoef between S & Cp and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));
    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);         % For the regression line
    x = [getIndexFns{fnn}(peakAlpha)(1) max(getIndexFns{fnn}(peakAlpha))];
    y = m*x + b;
     % Define bestIndexLabel here
    bestIndexLabel = sprintf("%s_{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
     % Print the slope, intercept, and standard error for each pair
    fprintf('Regression Results for %s and %s:\n', expData{2,edn}, bestIndexLabel);
    fprintf('Slope (m): %.4f, Intercept (b): %.4f\n', m, b);
    % Scatter plot
    this_figure = figure(3*(edn-1)+fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s^{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});
    % To find the standard error of fit
    % Calculate residuals
    residuals = expData{1, edn} - (m * getIndexFns{fnn}(peakAlpha) + b);
    % Calculate the sum of squared residuals (SSE)
    sse = sum(residuals.^2);
    % Calculate degrees of freedom (df)
    n = length(expData{1, edn});
    df = n - 2;  % Two parameters estimated: slope and intercept
    % Calculate the standard error of fit (SE)
    se = sqrt(sse / df);
    % Calculate the sum of squares of (x_i - mean(x)) for the x-values
    x_mean = mean(getIndexFns{fnn}(peakAlpha));
    x_variance = sum((getIndexFns{fnn}(peakAlpha) - x_mean).^2);
    % Standard error of the slope (SE_slope)
    se_slope = se / sqrt(x_variance);
    % Standard error of the intercept (SE_intercept)
    se_intercept = se * sqrt(1/n + (x_mean^2 / x_variance));
    fprintf('Standard Error of Fit for %s and %s: %.4f\n', ...
    bestIndexLabel, expData{2, edn}, se);
    fprintf('Slope (m): %.4f ± %.4f\n', m, se_slope);  % Standard error of slope
    fprintf('Intercept (b): %.4f ± %.4f\n', b, se_intercept);  % Standard error of intercept
    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");
    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)
    drawnow;
  end
end
if saveToFile
  % Save each figure to a separate file
  saveas(figure(1), "Scatter_ST1.jpg");
  saveas(figure(2), "Scatter_ST2.jpg");
  saveas(figure(3), "Scatter_CpT1.jpg");
  saveas(figure(4), "Scatter_CpT2.jpg");
end
