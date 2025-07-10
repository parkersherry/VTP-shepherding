posSeeds = zeros(15,1);
angle_seeds = zeros(15,1);

for i=1:15
load("/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Memory Vs Everything/mem_1_run_"+string(i)+".mat")
posSeeds(i) = (position_seed);
angle_seeds(i) = (angle_seed);
end

disp(posSeeds)

save("/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Memory Vs Everything/seeds.mat","angle_seeds","posSeeds")