% fname = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/FenceSuccess.mat';
fname = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/VaryVmax/results.mat';
load(fname)
% successRatio(successRatio>1) = 1;

fig1 = figure(1);
heatmap(nuArr,vmax,pols)
title("Phase Diagram: Nu and maximum sheep velocity vs order")
xlabel("Nu")
ylabel("Vmax")

% fig2 = figure(2);
% scatter(log2(alphaArr),successRatio,'filled')

