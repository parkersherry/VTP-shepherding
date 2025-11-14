dirname = "/Users/mikey/Summer 2025/MatlabGitVTP/MATLAB_code/SunOct5/";
nuArr = linspace(0.5,1,6);
memArr = linspace(1,120,120);
choice = 40;
% load(dirname+fname)
%
% blockMean = zeros(1,tmax);
% blockSize = memArr(choice);
%
% for t=1:tmax-blockSize
%     blockMean(t) = mean(Polarization_t(t:t+blockSize));
% end

% 
% 
tmin = 200;
tmaxB = 1000;
% BIG = zeros(3,120,tmaxB-tmin+1,6);
% for k=1:3
%     for nCounter =1:6
%     for i=1:120
%         fname = "nu_"+string(nuArr(6))+"_mem_"+string(memArr(i))+"_run_"+string(k)+".mat";
%         load(dirname+fname)
%         % X_Ti = X_T;
%         % fname = "nu_"+string(nuArr(6))+"_mem_"+string(memArr(i-1))+"_run_"+string(k)+".mat";
%         % load(dirname+fname)
%         BIG(k,i,:,nCounter) = ldod(tmin:tmaxB);
%     end
%     end
% end
% X = linspace(tmin,tmaxB,tmaxB-tmin+1);
% % toPlot = mean(BIG,1);
% figx = figure(125);
% scatter(X,squeeze(toPlot(:,1,:,1)),'filled','red')
% hold on
% scatter(X,squeeze(toPlot(:,end,:,1)),'filled','blue')

% fig1 = figure(1);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,1)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")
% 
% fig2 = figure(2);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,2)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")
% 
% fig3 = figure(3);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,3)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")
% 
% fig4 = figure(4);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,4)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")
% 
% fig5 = figure(5);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,5)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")
% 
% fig6 = figure(6);
% heatmap(X,memArr,squeeze(toPlot(:,:,:,6)))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% str = "$\Sigma_{i=1}^{N}\vert\vert X_T^m(t) - X_T^{m-1}\vert\vert $";
% title(" vs Polarization")


% fig2 = figure(2);
% toPlot(toPlot<100)=0;
% toPlot(toPlot>=100)=1;
% heatmap(linspace(1,1000,1000),memArr(1:end-1),squeeze(toPlot))
% xlabel("time (200 ms)")
% ylabel("Memory (s)")
% title("Thresholded difference in X_T")

% 
% N=35;
% 
% totalRuns = 15*numel(memArr)*numel(nuArr);
% counter = 0;
% 
% LDODS = zeros(120,6);
% pols = zeros(120,6);
% press = zeros(120,6);
% conv = zeros(120,6);
% for M=1:numel(memArr)
%     for U=1:numel(nuArr)
%         tempL = zeros(4,1);
%         tempPol = zeros(4,1);
%         tempPress = zeros(4,1);
%         tempConv = zeros(4,1);
%         for i=1:3
%             fname = "nu_"+string(nuArr(U))+"_mem_"+string(memArr(M))+"_run_"+string(i)+".mat";
%             load(dirname+fname)
%             tempPol(i) = mean(Polarization_t(200:end));
%             tempPress(i) = mean(pressures(200:end));
%             tempL(i) = mean(ldod(200:end));
%             tempConv(i) = mean(convexity_t(200:end));
%         end
%         LDODS(M,U) = mean(tempL);
%         pols(M,U) = mean(tempPol);
%         press(M,U) = mean(tempPress);
%         conv(M,U) = mean(tempConv);
%     end
% end

% bNum= 10;
% bSize = 120/bNum;
%
% convBins = zeros(bNum,6);
% pressBins = zeros(bNum,6);
% polsBins = zeros(bNum,6);
% ldodBins = zeros(bNum,6);
%
% for k=1:bNum
%     for nuT =1:6
%         convBins(k,nuT) = std(conv((k-1)*bSize+1:k*bSize,nuT));
%         pressBins(k,nuT) = std(press((k-1)*bSize+1:k*bSize,nuT));
%         polsBins(k,nuT) = std(pols((k-1)*bSize+1:k*bSize,nuT));
%         ldodBins(k,nuT) = std(LDODS((k-1)*bSize+1:k*bSize,nuT));
%     end
% end

% figConv = figure(1);
% hold on
%
% for j=1:6
%
%     scatter(linspace(1,bNum,bNum),convBins(:,j),'filled')
% end
% title("Convexity vs Memory")
% xlabel("Memory")
% ylabel("Convexity")
%
%
%
% figPols = figure(2);
% hold on
%
% for j=1:6
%
%     scatter(linspace(1,bNum,bNum),polsBins(:,j),'filled')
% end
% title("Polarization vs Memory")
% xlabel("Memory")
% ylabel("Polarization")
%
%
% figPress = figure(3);
% hold on
% for j=1:6
%
%     scatter(linspace(1,bNum,bNum),pressBins(:,j),'filled')
% end
% title("Pressure vs Memory")
% xlabel("Memory")
% ylabel("Pressure")
%
% figL = figure(4);
% hold on
%
% for j=1:6
%
%     scatter(linspace(1,bNum,bNum),ldodBins(:,j),'filled')
% end
% title("LDOD vs Memory")
% xlabel("Memory")
% ylabel("LDOD")

%%% SCATTER PLOTS
memMax = 120;
s = 6;
figConv = figure(1);
scatter(memArr(1:memMax),conv(1:memMax,s),'filled')
title("Convexity vs Memory")
xlabel("Memory")
ylabel("Convexity")



figPols = figure(2);
scatter(memArr(1:memMax),pols(1:memMax,s),'filled')

title("Polarization vs Memory")
xlabel("Memory")
ylabel("Polarization")


figPress = figure(3);
scatter(memArr(1:memMax),press(1:memMax,s),'filled')

title("Pressure vs Memory")
xlabel("Memory")
ylabel("Pressure")

figL = figure(4);
scatter(memArr(1:memMax),LDODS(1:memMax,s),'filled')
% 
title("LDOD vs Memory")
xlabel("Memory")
ylabel("LDOD")
% %
% %
% %
% % %HEAT MAPS
% %
%

% figConvheat = figure(5);
% heatmap(memArr,nuArr,conv')
% 
% title("Convexity Phase Diagram")
% xlabel("Memory")
% ylabel("Nu")
% 
% figPolsheat = figure(6);
% heatmap(memArr,nuArr,pols')
% 
% title("Polarization Phase Diagram")
% xlabel("Memory")
% ylabel("Nu")
% 
% 
% figPressheat = figure(7);
% heatmap(memArr,nuArr,press')
% 
% title("Pressure Phase Diagram")
% xlabel("Memory")
% ylabel("Nu")
% 
% figLheat = figure(8);
% heatmap(memArr,nuArr,LDODS')
% 
% title("LDOD Phase Diagram")
% xlabel("Memory")
% ylabel("Nu")
% % %
