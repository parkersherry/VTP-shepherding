tar = Target();
vmax = linspace(1/100,1/2,25);
tmaxarr = 5000;
% pols = cell(10,1);
% temp = zeros(10,1000);
nuArr = linspace(0,5,21);
temp = zeros(numel(vmax),numel(nuArr),5);

for seedIndex = 1:5
    for vmaxIndex = 1:25
        for nuIndex = 1:21



            fname = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/VaryVmax/vMax" + string(vmax(vmaxIndex)) +"_nu_"+string(nuArr(nuIndex))+ "_run_"+string(seedIndex)+".mat";
            load(fname)
            temp(vmaxIndex,nuIndex,seedIndex) = mean(Polarization_t(100:end));
        end
    end
end
pols = mean(temp,3);
disp(size(pols))
save("/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/VaryVmax/results.mat",'pols','nuArr',"vmax")