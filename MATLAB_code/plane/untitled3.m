Ndogs = 1;
memArr = 1:120;
Nsheep = 70;
fpath = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/lowMem/";
dtar = Target();
counter=0;
nuArr = linspace(0.5,1,10);


pols = zeros(5);
press = zeros(5);

for i=1:numel(nuArr)


    for j=1:numel(memArr)
        positionSeed = randi(100000);
        angSeed = randi(100000);

        seedPress = zeros(15,1);
        seedPols = zeros(15,1);
        for k=1:15
            fname = "nu_"+string(nuArr(i))+"_mem_"+string(memArr(j))+"_run_"+string(k)+".mat";
            load(fpath+fname)
            % scatter(1:numel(Polarization_t),Polarization_t)
            % counter=counter+1;
            tempPol = Polarization_t(200:end);
            tempPress = pressures(200:end);
            seedPress(k) = std(tempPress);
            seedPols(k) = std(tempPol);
        end
        pols(i,j) = mean(seedPols,'all');
        press(i,j) = mean(seedPress,'all');
    end
end

% for i=1:numel(press(:,1))
%     fig = figure(i);
%     scatter(1:numel(pols(i,:)),pols(i,:))
%     title("pressures, nu = "+string(nuArr(i)))
% end
% 
fig = figure(1);
heatmap(memArr,nuArr,pols)
title("Polarization Standard Deviation in Each Run")
hii = figure(2);
heatmap(memArr,nuArr,press)
title("Pressure Standard Deviation in Each Run")
