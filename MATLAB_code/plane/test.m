% %
% % N = ceil(linspace(50,250,12));
% % memory = ceil(linspace(1,240,12));
% %
% %
% % fenceTar = Target([20 20]);
% % quarterBubbleTar = Target([20 20]);
% % halfBubbleTar = Target([20 0]);
% % for memIndex = 1:numel(memory)
% %     disp("-------memIndex = " + string(memIndex) + '--------')
% %     for NIndex = 1:numel(N)
% %         disp("-------NIndex = " + string(NIndex) + '--------')
% %         for seedIndex = 1:5
% %
% %         posSeed = randi(100000);
% %         angSeed = randi(100000);
% %         fnameFenceHole = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune25/Fence/Fence_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
% %         fnameHalfBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune25/Half Bubble/HalfBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
% %         fnameQuarterBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune25/Quarter Bubble/QuarterBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
% %
% %
% %         serialSandbox(fnameFenceHole,N(NIndex),3.3,Inf,1000,fenceTar, ...
% %             memory(memIndex),posSeed,angSeed,1,'fence')
% %         serialSandbox(fnameHalfBubble,N(NIndex),3.3,Inf,1000,halfBubbleTar, ...
% %             memory(memIndex),posSeed,angSeed,1,'infiniteFence')
% %         serialSandbox(fnameQuarterBubble,N(NIndex),3.3,Inf,1000, ...
% %             quarterBubbleTar,memory(memIndex),posSeed,angSeed,1,'fenceNoGap')
% %         end
% %     end
% % end
%
%
%
% % N = 50;
% % memory = 240;
% % alpha = [10 100 1000 10000 10^5 10^6 10^7];
% %
% % fenceTar = Target([20 20]);
% % % quarterBubbleTar = Target([20 20]);
% % % halfBubbleTar = Target([20 0]);
% % for alphaIndex = 1:numel(alpha)
% %     disp("-------alphaIndex = " + string(alphaIndex) + '--------')
% %     for seedIndex = 1:5
% %
% %         posSeed = randi(100000);
% %         angSeed = randi(100000);
% %         fnameFenceHole = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Varying Alpha in Fence Situation High Speed/Fence_N_150_mem_240_alpha_'+string(alpha(alphaIndex)) + '_run_' + string(seedIndex) + '.mat';
% %         % fnameHalfBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune25/Half Bubble/HalfBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
% %         % fnameQuarterBubble = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/OvernightRunsJune25/Quarter Bubble/QuarterBubble_N_'+string(N(NIndex))+'_mem_'+string(memory(memIndex)) + '_run_' + string(seedIndex) + '.mat';
% %
% %
% %         serialSandbox(fnameFenceHole,N,3.3,alpha(alphaIndex),1000,fenceTar, ...
% %             memory,posSeed,angSeed,1,'fence')
% %         % serialSandbox(fnameHalfBubble,N(NIndex),3.3,Inf,1000,halfBubbleTar, ...
% %         %     memory(memIndex),posSeed,angSeed,1,'infiniteFence')
% %         % serialSandbox(fnameQuarterBubble,N(NIndex),3.3,Inf,1000, ...
% %         %     quarterBubbleTar,memory(memIndex),posSeed,angSeed,1,'fenceNoGap')
% %     end
% %
% % end
%
% load("/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Memory Vs EverythingN70/seeds.mat")
% N = [170 270];%ceil(linspace(25,300,20));
% memory = ceil(linspace(1,200,200));
% Nseeds = 15;
% tar = Target();
% for NIndex = 1:numel(N)
%     disp("---------------- " + newline+ newline+ newline  )
%     disp("NIndex = " + string(NIndex))
%     disp("---------------- "+ newline+ newline+ newline)
%     temp=zeros(numel(memory),Nseeds);
%     for seedIndex = 6:Nseeds
%         disp("-------seedIndex = " + string(seedIndex) + '--------')
%         posSeed = posSeeds(seedIndex);
%         angSeed = angle_seeds(seedIndex);
%         for memIndex = 1:numel(memory)
%             if(seedIndex == 6 && memIndex <56 && NIndex ==1)
%                 continue
%             end
%
%
%             fnameDog = '/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/Memory Vs EverythingN'+string(N(NIndex))+'/mem_' +string(memory(memIndex))+ '_run_'+string(seedIndex)+'.mat';
%
%
%             serialSandbox(fnameDog,N(NIndex),3.3,Inf,1500,tar,memory(memIndex),posSeed,angSeed,1,'zero');
%
%         end
%
%     end
% end
%
% aMean = mean(temp,2);
% aMin = min(temp,[],2);
% aMax = max(temp,[],2);
% xArr = linspace(1,15);
% modfun = @(a,x) a(1).*x+a(2);
% mdl = fitnlm(memory,aMean,modfun,[1/15 3]);
% aEstimated = mdl.Coefficients.Estimate;
% modFunT = @(t)modfun(aEstimated,t);
% regressPlt = arrayfun(modFunT,xArr);
% disp(mdl)
%
% fig1 = figure(1);
% hold on
% scatter(memory,aMean,'filled','black')
% scatter(memory,aMin,'filled','red')
% scatter(memory,aMax,'filled','blue')
% plot(xArr,regressPlt,'color','black','LineWidth',2)
% title("Timed Mixing vs Memory")
% save("/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/memVsMixExpDecay1_15.mat",'aMean','aMin','aMax')
% fig2 = figure(2);
% hold on
% scatter(memory,aMean,'filled','black')
% plot(xArr,regressPlt,'color','black','LineWidth',2)
% title("Timed Mean Mixing vs Memory")

tar = Target();
vmax = linspace(1/100,1/2,25);
tmaxarr = 5000;
% pols = cell(10,1);
% temp = zeros(10,1000);
nuArr = linspace(0,5,21);
disp(nuArr)
for seedIndex = 1:5
    disp(seedIndex)
    posSeed = randi(10000);
    angSeed = randi(10000);
    for vmaxIndex = 1:25
        for nuIndex = 1:21



            fname = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/VaryVmax/vMax" + string(vmax(vmaxIndex)) +"_nu_"+string(nuArr(nuIndex))+ "_run_"+string(seedIndex)+".mat";
            % load(fname)
            % temp(seedIndex,:) = Polarization_t;
            serialSandbox(fname,nuArr(nuIndex),50,3.3,Inf,tmaxarr,tar,1,posSeed,angSeed,0,'zero',vmax(vmaxIndex));
        end
    end
    % pols{vmaxIndex} = mean(temp,1);
end
% fig1 = figure(1);
% hold on
% for i=2:10
%     dname = 'vmax = ' + string(vmax(i));
%     plot(1:tmaxarr(i),pols{i},'DisplayName',dname,'LineWidth',i/2)
% end
% lgd = legend('FontSize',15,'Location','northeastoutside');
%
%
