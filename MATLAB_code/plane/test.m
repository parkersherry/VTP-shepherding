
N = ceil(linspace(50,250,22));
memory = ceil(linspace(1,240,22));
memory(2:end-1) = memory(2:end-1)-1;

fenceTar = Target([20 20]);
quarterBubbleTar = Target([20 20]);
halfBubbleTar = Target([20 0]);
for memIndex = 1:numel(memory)
    disp("-------memIndex = " + string(memIndex) + '--------')
    for NIndex = 1:numel(N)
        disp("-------NIndex = " + string(NIndex) + '--------')
        for seedIndex = 1:5

        posSeed = randi(100000);
        angSeed = randi(100000);
        fnameFenceHole = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune23/Fence/Fence_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
        fnameHalfBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune23/Half Bubble/HalfBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
        fnameQuarterBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune23/Quarter Bubble/QuarterBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
        

        serialSandbox(fnameFenceHole,N(NIndex),3.3,Inf,1000,fenceTar, ...
            memory(memIndex),posSeed,angSeed,1,'fence')
        serialSandbox(fnameHalfBubble,N(NIndex),3.3,Inf,1000,halfBubbleTar, ...
            memory(memIndex),posSeed,angSeed,1,'infiniteFence')
        serialSandbox(fnameQuarterBubble,N(NIndex),3.3,Inf,1000, ...
            quarterBubbleTar,memory(memIndex),posSeed,angSeed,1,'fenceNoGap')
        end
    end
end

