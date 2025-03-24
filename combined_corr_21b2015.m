% combined_corr_21b2015.m
% In this script, six values are closely approximated via golden section search
% - α value for which correlation coefficient ρ is strongest between S and T_1^a
% - α value for which correlation coefficient ρ is strongest between S and T_2^a
% - α value for which correlation coefficient ρ is strongest between Cp and T_1^a
% - α value for which correlation coefficient ρ is strongest between Cp and T_2^a
% Cp = Heat Capacity, S= Entropy
% T_1^a = First General Temperature Indices
% T_2^a = Second General Temperature Indices
% Additionally, curves ρ against α near these values are plotted in 2 figs
close all;     % Close any figures already opened
clear;         % and clear all variables
format long;   % More significant figures printed in the console
pkg load statistics;
% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = false; % Set to true to auto-save plots,
% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));
% Cell containing Entropy and Heat Capacity of 30 lower benzenoids hydrocarbon
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
   6 6 6 7 6  8 7  8 9 6 6  7  8  8  7  10 9  9  9  8  9  8  8 8 7  9  7  6  6  6  %(2,2)
   0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8  10 10 10 12 10 12 8 8 10 14 10 20 12 16 %(2,3)
   0 1 2 3 3  5 4  5 6 5 4  5  6  6  5  8  7  7  7  6  7  6  8 8 7  8  10 5  12 19 %(3,3)
 ]'; % Used for computing indices based on edge endpoint degree partitions
% Define v_exp here before using it in getIndexFns 30 lower BHs (v_exp: number of vertices)
v_exp = [
    6 10 14 14 18 18 18 18 18 16 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 26 22 26 24 32
]';
% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(alpha, v_exp) | n=1:T_1^{a} , n=2:T_2^{a}, a = alpha
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (v_exp - 2)), ((5 .* v_exp - 12) ./ ((v_exp - 2) .* (v_exp - 3))), (6 ./ (v_exp - 3))].^alpha, 2); % First General Temperature Indices (T_1^{a})
    @(alpha) sum(d_f .* [(4 ./ ((v_exp - 2).^2)), (6 ./ ((v_exp - 2) .* (v_exp - 3))), (9 ./ ((v_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_2^{a})
}';
% Cell containing their respective labels
indexName = {"T_1", "T_2"};
% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two S and Cp
numIndices = size(getIndexFns,2); % two (T_{1}^a  & T_{2}^a)
numCurves = numData*numIndices;   % two (4 curves)
% All x in visible ranges (both plots - near and far)
xnear = [linspace(-0.05, 0.21, 800); %for S
linspace(-0.04, 0.08, 800)]; %for Cp
xfar = linspace(-10,10,800); % for both S & Cp
% Do the same procedure for each experimental data i.e., Cp
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;
  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  %          ... These are coordinates where ρ-α of S-T_1^{a}, and Cp-T_1^{a}
  %          ... intersects ρ-α of S-T_2^{a} and Cp-T_2^{a}
   xmeet1 = [
    0.00080645;    % for S
    -0.00010702 % for Cp
  ](ii); % of these 2, either the first value or second is used, depending on ii
  ymeet1 = [
    0.91827; % for S
    0.97499  % for Cp
  ](ii);
  xmeet2 = [
    0.056853;    % for S
    0.0046624 % for Cp
  ](ii);
  ymeet2 = [
    0.91888; % for S
    0.975  % for Cp
  ](ii);
  ybox = [
  0.91895;
  0.97505
  % for figure 2 (blue dotted line)
  ](ii);
  % Plot the blue dashed box (before drawing the curves so it appear beneath)
  plot([xmeet1 xmeet1 xmeet2 xmeet2], [0 ybox ybox 0], '--b', 'LineWidth', lineWidth);
  yend = 0; % <-- to be assigned some value later for adjusting visible range
  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between S & Cp with specified α
    %                                and T_1^a /T_2^a (depending on n)
    get_indices_vals = @(alpha) getIndexFns{n}(alpha)(!isnan(expData{1,ii}));
    ccFn = @(alpha) corrcoef(
      get_indices_vals(alpha), % Either T_1  or T_2
      expData{1,ii}(!isnan(expData{1,ii})) % S or Cp
    )(1,2);
    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear(ii,:));
    yfar = arrayfun(ccFn, xfar);
    % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(GoldenSectionSearch_Maximum(ccFn, xnear(1), xnear(end), 1e-15));
    peakCorrCoeff = ccFn(peakAlpha);
    % Display peak alpha and peak correlation coefficient in the command window
    disp(['Peak Alpha: ', num2str(peakAlpha)]);
    disp(['Peak Correlation Coefficient: ', num2str(peakCorrCoeff)]);
    % Generate curve label S & Cp [T_1/T_2]
    curveLabels{n} = sprintf("%s and %s", expData{2,ii}, indexName{n});
    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;
    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear(ii,:), ynear, '-', 'LineWidth', lineWidth);
    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(ii,1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)
     yend = max(yend, ynear(end)); % y value to be used as visible y lower bound
  end
  % Mark and write on the plot the limits of alpha where T_1^alpha  is better than T_2^alpha
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark with blue asterisks
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]); % Write blue text
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);
  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "southeast"); % the location of the legend box
  ybox_space = [0.000005; % spacing between upper y-value with the blue dotted lines (figure 3)
                0.0005](ii);
  axis([xnear(ii,1) xnear(ii,end) yend ybox+ybox_space]); % adjust figure's visible range
  drawnow;
  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  set(leg, 'location', "southeast");
  if ii==2
    set(leg, 'location', "southeast");
  end
  hold off;
end
for ii = 1:4
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));
  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end
if saveToFile
  saveas(figure(1), "comb_corr_far_S.jpg");
  saveas(figure(2), "comb_corr_far_Cp.jpg");
  saveas(figure(3), "comb_corr_near_S.jpg");
  saveas(figure(4), "comb_corr_near_Cp.jpg");
end
