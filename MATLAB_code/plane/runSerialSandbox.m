dirname = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/SunOct5/";
N=100;
nuArr = linspace(0.5,1,6);
memArr = linspace(1,120,120);
totalRuns = 15*numel(memArr)*numel(nuArr);
counter = 0;
for i=1:5
    posSeed = randi(100000);
    angSeed = randi(100000);
    for M=1:numel(memArr)
        for U=1:numel(nuArr)

            fname = "nu_"+string(nuArr(U))+"_mem_"+string(memArr(M))+"_run_"+string(i)+".mat";
            serialSandbox(dirname+fname,nuArr(U),N,Inf,1000,Target([20,20]),memArr(M),posSeed,angSeed,1);
            disp("M = " + string(M) + "    |    i = " + string(i) + "    |    " + string(100*round(counter/totalRuns,4)) + "% completed")
            disp("-----------------------------------------------")
            counter=counter+1;
            

        end
    end
end