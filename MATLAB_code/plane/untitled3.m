clear all
% fname = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/FenceSuccess.mat';
fname = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/FenceSuccessAlphaw.mat';
load(fname)
successRatio(successRatio>1) = 1;

% fig1 = figure(1);
% heatmap(memory,NumAgents,successRatio)
% title("Phase Diagram: Memory and Number of Agents vs Success Ratio of Dog")
% xlabel("Memory Duration (# of time steps)")
% ylabel("Number of Agents")

fig2 = figure(2);
scatter(log2(alphaArr),successRatio,'filled')

