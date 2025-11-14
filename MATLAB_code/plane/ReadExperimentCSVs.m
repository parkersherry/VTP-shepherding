dirname = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/ExperimentalPosCSVs/Run";
Pressure = cell(16,1);
domainOfDanger = cell(16,1);
convexity = cell(16,1);
tmaxArr = [271,302,485,2029,1625,1901,1164,305,397,1071,90,550,212,306,274,305];

for runNum = 1:16
    disp(runNum)
    d = dirname+string(runNum)+"/";
    P = zeros(tmaxArr(runNum)-1,1);
    l = zeros(tmaxArr(runNum)-1,1);
    conv = zeros(tmaxArr(runNum)-1,1);
    for t=2:tmaxArr(runNum)
        fname = "output_t"+string(t)+".csv";
        T = readmatrix(d+fname,'HeaderLines', 1);
        agentIDs = T(:,1);
        N = numel(agentIDs);
        X = [T(:,1) T(:,2)];
        DT = delaunayTriangulation(X);
        P(t-1) = voronoiPressure(DT);
        l(t-1) = LDOD(X,75);
        [~,areaC] = convhull(X);
        shp = alphaShape(X);
        alpha = criticalAlpha(shp,'one-region');
        shp.Alpha = alpha;
        areaOfAlpha = area(shp);
        conv(t-1) = areaOfAlpha/areaC;
    end
    Pressure{runNum,1} = P;
    domainOfDanger{runNum,1} = l;
    convexity{runNum,1} = conv;
end

save("ExperimentalObservables.mat","convexity","domainOfDanger","Pressure")