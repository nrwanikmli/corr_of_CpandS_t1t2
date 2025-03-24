% good_alpha_21b2015.m
% In this script, 4 plots are generated which shows α interals for
% - good correlation coefficient ρ between S and T_1^a
% - good correlation coefficient ρ between S and T_2^a
% - good correlation coefficient ρ between Cp and T_1^a
% - good correlation coefficient ρ between Cp and T_2^a
% The limits of the intervals are also printed to console
close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 3;
fontSize = 16;
saveToFile = false; % Set to true to auto - save plots
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection, at the end of this script
% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));
% Cell containing Heat Capacity and Entropy of 30 lower benzenoids
expData = {reshape([% Entropy
   269.722 334.155 389.475 395.882 444.724 447.437 457.958 455.839 450.418 399.491
   499.831 513.857 508.537 507.395 506.076 512.523 500.734 520.307 509.210 513.879
   511.770 509.611 462.545 463.738 468.712 555.409 472.295 554.784 468.796 551.708
   ]', 30, 1), reshape([ % Heat Capacity
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
% Define nu_exp here before using it in getIndexFns 30 lower BHs (v_exp: number of vertices)
v_exp = [
    6 10 14 14 18 18 18 18 18 16 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 26 22 26 24 32
]';
% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(alpha, v_exp) | n=1:T_1^a, n=2:T_2^a
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (v_exp - 2)), ((5 .* v_exp - 12) ./ ((v_exp - 2) .* (v_exp - 3))), (6 ./ (v_exp - 3))].^alpha, 2); % First General Temperature Indices (T_1^a)
    @(alpha) sum(d_f .* [(4 ./ ((v_exp - 2).^2)), (6 ./ ((v_exp - 2) .* (v_exp - 3))), (9 ./ ((v_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_2^a)
}';
% Cell containing their respective labels
indexName = {"T_1", "T_2"};
% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two S and Cp
numIndices = size(getIndexFns,2); % two (T_1  & T_2)
numCurves = numData*numIndices;   % two (4 curves)
% Boundaries for visible intervals, for each index-property pair
%             T_1     T_2
xstart = [   -0.5     -0.3       % S
             -0.5     -0.2];     % Cp
xend = [      0.5      0.25      % S
              0.3      0.1];     % Cp
ystart = [    0.92455  0.925     % S
              0.9773   0.9775];  % Cp
yend = [      0.91     0.91      % S
              0.972    0.973];   % Cp
% Exact rho value considered good for
%              S       Cp
a_goodrho = [0.91885, 0.9750];
% Colors (different shades of cyan and green)
colShaded = {[0.85, 1, 1]; [0.85, 1, 0.85]};
colIndicatorVert = {[0.2, 0.55, 0.55]; [0, 0.5, 0]};
colIndicatorHorz = {[0.35, 0.75, 0.75]; [0.45, 0.7, 0.45]};
colCurve = {[0, 0.75, 0.75]; [0, 0.5, 0]};
  % Do the same procedure for each experimental data i.e., Cp, S
for ii = 1:numData
   % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between S & Cp with specified α
    %                          and T_{1}^a & T_{2}^a (depending on n)
    ccFn = @(alpha) corrcoef(
    getIndexFns{n}(alpha)(!isnan(expData{1,ii})), % Either T_1 or T_2
    expData{1,ii}(!isnan(expData{1,ii})) % S or Cp
    )(1,2);
     this_fig = figure('Name', sprintf('Correlation between %s and %s', indexName{n}, expData{2,ii})); % Adjusted line
     hold on;
    % Get Interval limits and print them to console
    % Get alpha for highest rho first
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -4, 0, 1e-15));
    % Calculate correlation coefficient (rho) for the peakAlpha
    rho_value = ccFn(peakAlpha);
    % Prepare a funcion to calc |rho(a)-goodrho|
    ccFn_good = @(alpha)(-abs(ccFn(alpha)-a_goodrho(ii)));
    % and func to get the limit, i.e., value of alpha where rho is 0.91885 or 0.9750
    getLimitFromInterval = @(lb, ub) mean(
      GoldenSectionSearch_Maximum(ccFn_good, lb, ub, 1e-15));
    a_lb = getLimitFromInterval(peakAlpha-1, peakAlpha); % Search to the left
    a_ub = getLimitFromInterval(peakAlpha, peakAlpha+1); % Search to the right
    % Write the intervals in console
    disp(sprintf("ρ(%s,%s_α) ≥ %.02f when α ∈ [%.08f, %.08f]",
         expData{2,ii}, indexName{n}, a_goodrho(ii), a_lb, a_ub));
    % Display the rho value for the current alpha
    fprintf('For %s and %s: α = %.4f, ρ = %.6f\n', expData{2,ii}, indexName{n}, peakAlpha, rho_value);
    % Plot the actual curve, but exclude good alpha range (<- separately drawn)
    % generate x values
    x = [linspace(xstart(ii,n),a_lb,500), linspace(a_ub,xend(ii,n),500)];
    % generate the corresponding y values
    y = arrayfun(ccFn, x);
    plot(x, y, '-', 'LineWidth', lineWidth);
    % Shade the area inside the good alpha interval
    u0 = a_lb;         u_width = a_ub-a_lb;
    v0 = ystart(ii,n); v_height = yend(ii,n) - ystart(ii,n);
    rectangle('Position', [u0, v0, u_width, v_height], 'FaceColor', colShaded{ii}, 'LineStyle', 'none');
    % Draw the indicator lines for the good alpha interval
    % Vertical dashed lines:
    plot([a_lb a_lb], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([a_ub a_ub], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    % Horizontal black dashed line, horizontal colored dashed line:
    plot([xstart(ii,n) a_lb], [a_goodrho(ii) a_goodrho(ii)],
         '--k', 'LineWidth', lineWidth/1.75);
    plot([a_lb a_ub], [a_goodrho(ii) a_goodrho(ii)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});
    % Write on the plot the interval limits
    text(a_lb, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_lb)))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    text(a_ub, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_ub)))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    % Finally, plot the colored curve within the good alpha range
    x_in = linspace(a_lb,a_ub,400);
    y_in = arrayfun(ccFn, x_in);
    plot(x_in, y_in, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});
    % Also highlight the good interval on the x-axis
    plot([a_lb, a_ub], [yend(ii,n), yend(ii,n)], '-',
         'LineWidth', lineWidth, 'Color', colCurve{ii});
    % Label the plot                                [S, Cp ]        [T^1,T^2]
    title(sprintf('Correlation between %s and %s', expData{2,ii}, indexName{n}));
    xlabel('α');
    ylabel('ρ');
    drawnow;
    % Now apply the axis limits with dynamic bounds
    axis([xstart(ii,n) xend(ii,n) yend(ii,n) ystart(ii,n)])
    % Replace all hypens with minuses
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    % Set the font size
    set(findall(this_fig,'-property','FontSize'),'FontSize', fontSize);
    drawnow;
    hold off;
  end
end
if saveToFile
  saveas(figure(1), "good_alpha_ST1.jpg");
  saveas(figure(2), "good_alpha_ST2.jpg");
  saveas(figure(3), "good_alpha_CpT1.jpg");
  saveas(figure(4), "good_alpha_CpT2.jpg");
end
