fpath = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/lowMem/";
pols = zeros(20,1);
press = zeros(20,1);
memArr = ceil(linspace(1,60,20));


for i=1:20
    seedPress = zeros(4,1);
    seedPols = zeros(4,1);
    for k=1:4

        fname = "_mem_"+string(memArr(i))+"_run_"+string(k)+".mat";
        load(fpath+fname)
        tempPol = Polarization_t(200:end);
        if (i==1 && k==1)
            figggg = figure(10);
            scatter(1:numel(tempPol),tempPol)
        end
        tempPress = pressures(200:end);
        seedPress(k) = mean(tempPress);
        seedPols(k) = mean(tempPol);

    end
    pols(i) = mean(seedPols,'all');
    press(i) = mean(seedPress,'all');
end

fig = figure(1);
scatter(memArr,pols,'filled')
title("polarization vs memory duration")

hii = figure(2);
scatter(memArr,press,'filled')
title("pressure vs memory duration")
