numruns=10050;
tmax = 1000;
tDomain = 1:tmax;

% nu = zeros(numruns,1);
% nustr = ["04","08","12","16"];
% Ms = zeros(tmax,numruns);
% ms = zeros(tmax,numruns);
% rmed = zeros(tmax,numruns);
% tri_root_area = zeros(tmax,numruns);
%CM = zeros(tmax,numruns,2);
% angmoms = zeros(tmax,numruns);
% absangmoms = zeros(tmax,numruns);

N = [50 75 80 105 110 135 140 165 170 195 200 225 230 255 260 285 290 315 320 345];
nu = [0.001 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5 9.75 10];

All_polarization = zeros(tmax,numruns);
Pressures = zeros(tmax,numruns);
Drift = zeros(tmax,numruns);

Xruns = zeros(tmax,numruns);
Yruns = zeros(tmax,numruns);
i=0;

IsoPolarizations = zeros(numel(nu),numel(N));
tempMeanPols = zeros(5,1);

IsoDrifts = zeros(numel(nu),numel(N));
tempMaxDrift = zeros(5,1);

mixing = zeros(numel(nu),numel(N));

Nstr = string(N);
nustr = string(nu);
for nuIndex=1:numel(nu)
    %disp(nuIndex)
    for NIndex=1:20
        if (mod(nuIndex,2)==mod(NIndex,2))
            temp = 0;
            for seedIndex=1:5
                i=i+1;
                filename = '/Users/mikey/Summer 2025/plane/OvernightM26/Dog/D_N_'+Nstr(NIndex) + '_nu_'+nustr(nuIndex) + '_run_' + string(seedIndex) + '.mat';
                load(filename);
                X_0 = DT_t{1}.Points;
                X_f = DT_t{tmax}.Points;
                DT_forMix = delaunayTriangulation(X_f(Ndogs+1:end,:));
                temp = temp + unifMixMetric(N,Ndogs,DT_forMix,X_0);
                %pol_t = zeros(tmax,1);
                %angmom_t = zeros(tmax,1);
                %absangmom_t = zeros(tmax,1);
                %A_t = zeros(tmax,1);


                % for t=1:tmax
                %     DT = DT_t{t};
                %     X = DT.Points;
                %     %U = U_t(:,:,t);
                %     CM = mean((X(Ndogs+1:end,:)));
                %     Drift(t,i) = vecnorm(CM);
                %     % Xruns(i,t) = CM(1);
                %     % Yruns(i,t) = CM(2);
                %     %CM(t,i,:) = mean(X(2:end,:));
                %     %DogDistFromCM(t,i) = vecnorm(CM(t,i,:) - X(1,:),2,2);
                %     %pol_t(t) = polarization(U,Ndogs);
                %     %[angmom_t(t),absangmom_t(t)] = angularMomentum(X,U);
                %     %A_t(t) = inwardTotalArea(DT);
                %     %com = mean(X);
                %     %rmed(t,i) = median(vecnorm(X-com,2,2));
                %     %tri = pointLocation(DT,com);
                %     %verts = DT.ConnectivityList(tri,:);
                %     %coords = X(verts,:);
                %     %tri_root_area(t,i) = sqrt(poly_area(coords(:,1),coords(:,2)));
                %     %Pressures(t,i) = voronoiPressure(DT);
                %     All_polarization(t,i) = polarization(U_t(:,:,t),Ndogs);
                % end
                % tempMeanPols(seedIndex) = mean(All_polarization(201:end,i));
                % tempMaxDrift(seedIndex) = max(Drift(:,i));

            end
            temp = temp/5;
            mixing(nuIndex,NIndex) = temp;
            % IsoPolarizations(nuIndex,NIndex) = mean(tempMeanPols);
            % IsoDrifts(nuIndex,NIndex) = mean(tempMaxDrift);
            % disp([nuIndex NIndex])
        end
    end
    disp(nuIndex);
end

save('/Users/mikey/Summer 2025/plane/SheepObservedPolarizationAndDrifts1.mat',"mixing",'nu','N');
% MeanPressure = mean(Pressures,2);
% %rmm = mean(rmed,2);
% plotDist = mean(DistMax,2);
%
% fig2 = figure(Name="Mean Pressure Graph");
% scatter(tDomain, MeanPressure,'filled');
% xlabel('Number of Time Steps (Arb. Units)')
% ylabel('Mean Voronoi Pressure [1/area]')
% title("Mean Pressure vs Time")
% legend("yaaass",'Location','northeastoutside')


% fig5 = figure(name="Trajectory of CM");
% x = tDomain;
% y = sin(tDomain);
% plot(Xruns(1,:), Yruns(1,:)), grid on;
% xlabel('X');
% ylabel('Y');
% title('Trajectory of a Point');
