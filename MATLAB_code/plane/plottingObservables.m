filename = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/plane/SheepObservedPolarizationAndDrifts.mat';
load(filename);

clear("fig1")
clear("fig2")
clear("fig3")
clear("fig4")
clear("fig5")
clear("fig6")
clear("fig7")
clear("fig8")
clear("fig9")


N = [50 75 80 105 110 135 140 165 170 195 200 225 230 255 260 285 290 315 320 345];
nu = [0.001 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5 9.75 10];

oddsN = linspace(1,19,10);
evensN = linspace(2,20,10);

oddsNu = linspace(1,41,21);
evensNu = linspace(2,40,20);


%---------------------------------%

fig1 = figure(1);
hold on

scatter(nu(evensNu),IsoDrifts(evensNu,evensN),"filled")
scatter(nu(oddsNu),IsoDrifts(oddsNu,oddsN),"filled")

hold off
xlabel('Alignment Strength \nu') 
ylabel('Mean Maximum Drift') 
legend('N = ' + string(N),'location','northeastoutside');

title("Maximum Drift vs N")
%---------------------------------%

fig2 = figure(2);
hold on

scatter(N(evensN),IsoDrifts(evensNu,evensN),"filled")
scatter(N(oddsN),IsoDrifts(oddsNu,oddsN),"filled")

hold off

legend('\nu = ' + string(nu),'location','northeastoutside');
xlabel('Number of Sheep') 
ylabel('Mean Maximum Drift') 
title("Maximum Drift vs N")

%---------------------------------%

fig3 = figure(3);
hold on

scatter(nu(evensNu),IsoPolarizations(evensNu,evensN),"filled")
scatter(nu(oddsNu),IsoPolarizations(oddsNu,oddsN),"filled")

hold off
xlabel('Alignment strength \nu') 
ylabel('Mean Polarization') 
legend('N = ' + string(N),'location','northeastoutside');

title("Mean Non-Transient Polarization vs \nu")
%---------------------------------%

fig4 = figure(4);
hold on

scatter(N(evensN),IsoPolarizations(evensNu,evensN),"filled")
scatter(N(oddsN),IsoPolarizations(oddsNu,oddsN),"filled")

hold off

legend('\nu = ' + string(nu),'location','northeastoutside');
title("Mean Non-Transient Polarization vs N")
xlabel('Number of Sheep') 
ylabel('Mean Polarization') 


filename1 = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/plane/SheepObservedPolarizationAndDrifts1.mat';
load(filename1);

N = [50 75 80 105 110 135 140 165 170 195 200 225 230 255 260 285 290 315 320 345];
nu = [0.001 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5 9.75 10];
%---------------------------------%
fig5 = figure(5);
hold on

scatter(N(evensN),mixing(evensNu,evensN),"filled")
scatter(N(oddsN),mixing(oddsNu,oddsN),"filled")
title("Mean Mixing vs N")
xlabel('Number of Sheep') 
ylabel('Mixing Metric') 
hold off
% 
legend('\nu = ' + string(nu),'location','northeastoutside');
title("Mean Mixing vs N")
%---------------------------------%
fig6 = figure(6);
hold on

scatter(nu(evensNu),mixing(evensNu,evensN),"filled")
scatter(nu(oddsNu),mixing(oddsNu,oddsN),"filled")

hold off

legend('N = ' + string(N),'location','northeastoutside');
xlabel('Alignment Strength \nu') 
ylabel('Mean Mixing Parameter') 
title("Mean Mixing vs \nu")
%---------------------------------%
fig7 = figure(7);
hold on

for k=1:numel(evensN)
    scatter(nu(evensNu),N(evensN(k)).*ones(size(nu(evensNu))),[],mixing(evensNu,evensN(k)),'filled')
end
for k=1:numel(oddsN)
    scatter(nu(oddsNu),N(oddsN(k)).*ones(size(nu(oddsNu))),[],mixing(oddsNu,oddsN(k)),'filled')
end
xlabel('Alignment Strength \nu') 
ylabel('Number of Sheep') 
title("Mean Mixing Phase Plot")
cb = colorbar();
xlim([-0.5 10.5])
ylim([30 350])
grid on

hold off

fig8 = figure(8);
hold on

for k=1:numel(evensN)
    scatter(nu(evensNu),N(evensN(k)).*ones(size(nu(evensNu))),[],IsoPolarizations(evensNu,evensN(k)),'filled')
end
for k=1:numel(oddsN)
    scatter(nu(oddsNu),N(oddsN(k)).*ones(size(nu(oddsNu))),[],IsoPolarizations(oddsNu,oddsN(k)),'filled')
end
xlabel('Alignment Strength \nu') 
ylabel('Number of Sheep') 
title("Mean Non-Transient Polarization Phase Plot")
cb = colorbar();
xlim([-0.5 10.5])
ylim([30 350])
grid on

hold off


fig9 = figure(9);
hold on

for k=1:numel(evensN)
    scatter(nu(evensNu),N(evensN(k)).*ones(size(nu(evensNu))),[],IsoDrifts(evensNu,evensN(k)),'filled')
end
for k=1:numel(oddsN)
    scatter(nu(oddsNu),N(oddsN(k)).*ones(size(nu(oddsNu))),[],IsoDrifts(oddsNu,oddsN(k)),'filled')
end
xlabel('Alignment Strength \nu') 
ylabel('Number of Sheep') 
title("Mean Maximum Drift Phase Plot")
cb = colorbar();
xlim([-0.5 10.5])
ylim([30 350])
grid on

hold off
