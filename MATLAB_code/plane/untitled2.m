Ndogs = 1;
memArr = 61:120;
Nsheep = 70;
fpath = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/lowMem/";
dtar = Target();
counter=0;
nuArr = linspace(0.5,1,10);
total=15*numel(memArr)*numel(nuArr);
for k=1:15
    load(fpath+"nu_0.5_mem_1_run_"+string(k)+".mat")
    for i=1:numel(memArr)
        for j=1:numel(nuArr)
            disp("")
            disp("-----------------------------------------------------------")
            str = "k = "+string(k)+"    |    i = " + string(i)+"    |    j = "+string(j) + "    |    "+string(round(10000*counter/total,2)/100) + "% completed";
            disp(str)
            fname = "nu_"+string(nuArr(j))+"_mem_"+string(memArr(i))+"_run_"+string(k)+".mat";
            a = serialSandbox(fpath+fname,nuArr(j),Nsheep,Inf,1000,dtar,memArr(i),position_seed,angle_seed,Ndogs);
            counter=counter+1;
        end

    end
end


