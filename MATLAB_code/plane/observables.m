clear all
N = 70;
memory = ceil(linspace(1,15,15));
Nseeds = 5;
times = linspace(1,1500,1500);

convexity = zeros(numel(memory),Nseeds);
pressure = zeros(numel(memory),Nseeds);
polarizationArr = zeros(numel(memory),Nseeds);
for seedIndex = 1:Nseeds
    disp("-------seedIndex = " + string(seedIndex) + '--------')
    for memIndex = 1:numel(memory)
        fname = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Memory Vs Everything/mem_' +string(memory(memIndex))+ '_run_'+string(seedIndex)+'.mat';
        load(fname);
        areasRatio = zeros(1,tmax);
        Polarization_t = zeros(1,tmax);
        pressures_t = zeros(1,tmax);
        for t=1:tmax 
            [~,ConvArea] = convhull(X_T(:,:,t));
            shp = alphaShape(X_T(:,:,t));
            alpha = criticalAlpha(shp,'one-region');
            shp.Alpha = alpha;
            areasRatio(t) = area(shp)/ConvArea;
            Polarization_t(t) = polarization(U_t(:,:,t),1);
            pressures_t(t) = voronoiPressure(delaunayTriangulation(X_T(:,:,t)));
        end
        convexity(memIndex,seedIndex) = mean(areasRatio(50:end));
        pressure(memIndex,seedIndex) = mean(pressures_t(50:end));
        polarizationArr(memIndex,seedIndex) = mean(Polarization_t(50:end));
    end
end

convexity = mean(convexity,2);
pressure = mean(pressure,2);
polarizationArr = mean(polarizationArr,2);

fig1 = figure(1);
scatter(memory,convexity,'filled','blue')

fig2 = figure(2);
scatter(memory,pressure,'filled','blue')

fig3 = figure(3);
scatter(memory,polarizationArr,'filled','blue')