load("ExperimentalObservables.mat")
figPress = figure(1);
hold on
tmaxArr = [271,302,485,2029,1625,1901,1164,305,397,1071,90,550,212,306,274,305];


for runNum = 1:3
    tm = tmaxArr(runNum);
    scatter(linspace(1,tm-1,tm-1),Pressure{runNum},'filled')
end