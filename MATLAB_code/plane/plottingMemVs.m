%clear all
load('/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/memVsAllN170.mat')


fig1 = figure(1);
scatter(memory,pressure,'blue','filled')
title("Memory vs Pressure")
xlabel("Memory Duration")
ylabel("Non-Transient Pressure")

fig2 = figure(2);
scatter(memory,polarizationArr,'blue','filled')
title("Memory vs Polarization")
xlabel("Memory Duration")
ylabel("Non-Transient Polarization")

fig3 = figure(3);
scatter(memory,convexity,'blue','filled')
title("Memory vs Convexity")
xlabel("Memory Duration")
ylabel("Non-Transient Convexity")


fig4 = figure(4);
scatter(memory(1:40),pressure(1:40),'blue','filled')
title("Memory vs Pressure in Low Memory Regime")
xlabel("Memory Duration")
ylabel("Non-Transient Pressure")

fig5 = figure(5);
scatter(memory(1:40),polarizationArr(1:40),'blue','filled')
title("Memory vs Polarization in Low Memory Regime")
xlabel("Memory Duration")
ylabel("Non-Transient Polarization")

fig6 = figure(6);
scatter(memory(1:40),convexity(1:40),'blue','filled')
title("Memory vs Convexity in Low Memory Regime")
xlabel("Memory Duration")
ylabel("Non-Transient Convexity")