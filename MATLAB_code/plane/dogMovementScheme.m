%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOGSWEEPVORONOI
% Computes and specifies the dynamical rules adhered to by dog agents.
%
% Syntax:
%   dSV = dogSweepVoronoi(X, DT, Ndogs, L, tar, firstNeighbors)
% Input:
%   X               :   (:,2), Positions of agents
%   DT              :   delaunayTriangulation object
%   Ndogs           :   Number of dogs
%   L               :   positive double.
%                       Effective length scale of repulsion.
%   tar             :   Target object. Used to find the equilibrium
%                       position of dog agents using findDogEquil.
%   firstNeighbours :   Nx1 cell array of all agents'
%                       first Voronoi neighbours. Use nbhd from
%                       neighbourhoods!
% Output:
%   P           :   Ndogs x 2 double. P(i,:) is the displacement vector
%                   for the ith dog agent.
%   indexMax    :   Ndogs x 1 integer. indexMax(i) is the index
%                   of the sheep agent which maximizes the angle
%                   between the vector from the ith dog agent to the local
%                   CM and the vector from the ith dog to the maximizing
%                   agent.
%   positions   :   Ndogs x 2 double. positions(i,:) is the position
%                   the ith dog agent wants to move to. positions(i,:)
%                   is defined such that the ith dog agent and its target
%                   sheep agent are collinear with the local CM.
function output = dogMovementScheme(X_T, U,DT, Ndogs, L, tar ,t,dt ,LastSeenCell,XmemCell,scalarField,prefVel,alpha,memDuration)

arguments
    X_T
    U
    DT
    Ndogs % int8 {mustBePositive}
    L double {mustBePositive}
    tar
    t
    dt
    LastSeenCell
    XmemCell
    scalarField
    prefVel
    alpha
    memDuration = 240
end
% toggles whether to prioritize being collinear with
% local CM and dog target
stopSimul = false;
advectRememberedAgents = false;
X = X_T(:,:,t);

N = numel(X)/2;
Xsheep = X(Ndogs+1:end,:);

recallDelay = memDuration;

[nbhd, nearest, ~] = neighborhoods(DT,2);

shp = alphaShape(Xsheep);
alpha = criticalAlpha(shp,'one-region');
shp.Alpha = alpha;
[~,P] = boundaryFacets(shp);
areaOfAlpha = area(shp);
if (~isempty(P))
    P(end+1,:) = P(1,:);
end
[indicesConvHull,areaConvHull] = convhull(Xsheep);
ConvHullTarget = Target(Xsheep(indicesConvHull,:));
CMs = zeros(Ndogs,2);

% keeps track of sheep agent indices while finding maximizing agent
SheepToCollectArr = zeros(Ndogs,1);

% initialize dog agent displacement vector
DogDisplacementVec = zeros(Ndogs, 2);

% initialize positions to which dog agents move
goalLocationArr = DogDisplacementVec;

% stalk = zeros(Ndogs,1);
% dogStalkLength = 3.5;
%%%%%%%%%%THIS ONLY WORKS FOR ONE DOG

exponentialDecay = zeros(N,2);
decayRate = recallDelay*0.02;
ExpAmplitude = 0.43;
ExpShift = 0.48;

repelConvHullStrength = 1;

sheepIsSeen = ones(Ndogs,1);

for i = 1:Ndogs


    if t*dt < recallDelay
        timeThreshold = 0;
    else
        timeThreshold = t*dt - recallDelay;
    end

    SheepNbhdPast = find(LastSeenCell{i} > timeThreshold);
    forgotten = find((~(LastSeenCell{i} > timeThreshold)).*(LastSeenCell{i}~=0));
    

    % for j=1:J
    %     X_j = X_T(:,:,LastSeen(SheepNbhdPast(j)));
    %     TotalX(SheepNbhdPast(j),:) = X_j(SheepNbhdPast(j),:);
    % end

    %-----------------------%
    %if you can see the residual image of the sheep then obvi there is no
    %sheep there and you should forget it
    if (~isempty(SheepNbhdPast))
        ToDel = zeros(size(SheepNbhdPast));
        inside = inhull(XmemCell{i}(SheepNbhdPast,:),P);
        for j=1:numel(SheepNbhdPast)
            if (inside(j)==1)
                continue
            end
            Xcoords = [X(i,1) XmemCell{i}(SheepNbhdPast(j),1)];
            %disp(size(Xcoords))
            Ycoords = [X(i,2) XmemCell{i}(SheepNbhdPast(j),2)];
            [xIntersect,~] = polyxpoly(Xcoords,Ycoords,P(:,1),P(:,2));
            if (isempty(xIntersect))
                ToDel(j) = 1;
            end
        end
        ToDel = find(ToDel);
        XmemCell{i}(SheepNbhdPast(ToDel),:) = 0;
        LastSeenCell{i}(SheepNbhdPast(ToDel)) = 0;
    
        SheepNbhdPast(ToDel) = [];
    end
    
    %-----------------------%
    CurrAllNbhd = [nbhd{i,1} nbhd{i,2}];

    CurrAllNbhd = unique(CurrAllNbhd);
    forgotten = setdiff(forgotten,CurrAllNbhd);
    CurrSheepnbhd = CurrAllNbhd(find(CurrAllNbhd>Ndogs));

    SubflockCells = {Xsheep};
    % SubflockCells = getSubflocks(Xsheep,X(CurrSheepnbhd,:),L,N,Ndogs,{},0);

    XmemCell{i}(CurrSheepnbhd,:) = X(CurrSheepnbhd,:);
    
    LastSeenCell{i}(CurrSheepnbhd) = t*dt;
    currMeanVel = mean(U(CurrSheepnbhd,:),1);

    %TotalX(CurrSheepnbhd,:) = X(CurrSheepnbhd,:);
    SheepNbhd = unique(union(SheepNbhdPast,CurrSheepnbhd));
    
    XmemCell{i}(forgotten,:) = 0;
    if strcmp(scalarField,'fence')
        indices = find(~((XmemCell{i}(SheepNbhd,1)<20).*(XmemCell{i}(SheepNbhd,2)<20)));
        SheepNbhd(indices) = [];
        if (isempty(SheepNbhd))
            disp("all sheep are outside the fence")
            stopSimul = true;
            break
        end
    end
    
    SheepNbhd( find(~any(XmemCell{i}(SheepNbhd,:)')), : ) = [];

    if (isempty(SheepNbhd))
        sheepIsSeen(i) = 0;
        continue
    end

    exponentialDecay = ExpAmplitude.*exp(decayRate.*(LastSeenCell{i}(SheepNbhd)-t*dt)./recallDelay) + ExpShift;
    sumExpDecay=(sum(exponentialDecay));
    CM = sum(XmemCell{i}(SheepNbhd,:).*exponentialDecay,1)./sumExpDecay;
    equilibrium = findDogEquil(CM,L , tar,numel(SheepNbhd));

    % finds sheep agent with indexMax

    xDog = X(i,:);
    CMtoDog = xDog - CM;
    CMtoDogNorm = vecnorm(CMtoDog,2,2);
    if (CMtoDogNorm~=0)
        CMtoDog = CMtoDog ./CMtoDogNorm ;
    end

    DogToSheep = xDog-XmemCell{i}(SheepNbhd,:);
    DogToSheepNorm = vecnorm(DogToSheep,2,2);
    if (DogToSheepNorm~=0)
        DogToSheep = DogToSheep./DogToSheepNorm;
    end

    thetaArray = acos(DogToSheep*transpose(CMtoDog));

    %tracks the sheep being targed by other dogs
    [~,~,intersectionIndices] = intersect(SheepToCollectArr,SheepNbhd);

    %ensures this dog doesn't target those sheep
    thetaArray(SheepNbhd(intersectionIndices)) = 0;
    [~,tempIndex] = max(thetaArray);%.*exponentialDecay);
    SheepToCollect = SheepNbhd(tempIndex);

    %ensures all other dogs don't target this sheep
    SheepToCollectArr(i) = SheepToCollect;


    %repel from concave hull
    h = homeToTarget(ConvHullTarget,X(i,:),1,false);
    distToConvHull = vecnorm(h,2,2);
    if (distToConvHull ~= 0)
        h = h./distToConvHull;
    end

    % equilibrium ON protoco


    XsheepMaxIndex = XmemCell{i}(SheepToCollect,:);

    goalLocation = (XsheepMaxIndex - CM);
    goalLocationNorm = vecnorm(goalLocation,2,2);
    if (goalLocationNorm~=0)
        goalLocation =  goalLocation./goalLocationNorm;
        goalLocation = goalLocation.*L;
    end
    goalLocation = goalLocation + XsheepMaxIndex;
    dogToGoal = goalLocation - xDog;

    % if vecnorm(dogToGoal,2,2) < dogStalkLength
    %     stalk(i) = 1;
    %     continue
    % end
    dogToGoalNorm = vecnorm(dogToGoal,2,2);
    if (dogToGoalNorm~=0)
        dogToGoal = dogToGoal ./dogToGoalNorm;
    end
    
    %---------------
    CMtoDogTar = homeToTarget(tar,CM);
    CMtoDogTarNorm = vecnorm(CMtoDogTar,2,2);
    if (CMtoDogTarNorm~=0)
        CMtoDogTar = CMtoDogTar./CMtoDogTarNorm;
    end
    thetaDogEquil = acos(dot(CMtoDog,CMtoDogTar));


    DogToEquil = equilibrium - xDog;
    DogToEquilNorm = vecnorm(DogToEquil,2,2);
    % if DogToEquilNorm < dogStalkLength
    %         stalk(i) = 1;
    %         continue
    % end
    if (DogToEquilNorm~=0)
        DogToEquil = DogToEquil./DogToEquilNorm;
    end

    s = transition(distToConvHull/L, 'expReciprocal');

    if tar.emptiness
        sDriving = 0;
    else
        sDriving = transition(thetaDogEquil/pi,'expReciprocal');
    end
    nearestIndex = nearest(i);
    dog_nearestRepel = (xDog - X(nearestIndex,:));

    dog_nearestRepelNorm = vecnorm(dog_nearestRepel,2,2);
    if (dog_nearestRepelNorm~=0)
        dog_nearestRepel = dog_nearestRepel./dog_nearestRepelNorm;
    end
    % dogS = 0;
    % 
    % if nearest(i)<Ndogs+1
    % 
    % end
    %turn off repel from CH
    if nearest(i)<Ndogs+1
        dogS = 0;
    else
        dogS = transition(dog_nearestRepelNorm/L, 'expReciprocal');
    end
    
    %

    DogDisplacementVec(i,:) =  (1-dogS) .* dogToGoal.*(1-sDriving)+ DogToEquil.*(sDriving);
    DogDisplacementVec(i,:) = DogDisplacementVec(i,:) + (dogS).*dog_nearestRepel-(0*s*repelConvHullStrength).* h;
    goalLocationArr(i,:) = goalLocation;

    CMs(i,:) = CM;
    %Xtotals{i} = Xmem;
    if advectRememberedAgents
        XmemCell{i}(SheepNbhd,:) = XmemCell{i}(SheepNbhd,:) + currMeanVel;
    end
end

if (strcmp(scalarField,'fence'))
    U(1:Ndogs,:) = U(1:Ndogs,:)+prefVel;
end


if (max(sheepIsSeen) == 0)
    disp("No Dog Sees Any Sheep -- Stopping Simulation")
    stopSimul = true;
end



clear("chullTar")
W = (DogDisplacementVec);%+U(1:Ndogs,:);
g = ones(Ndogs,1);
tooLarge = find(vecnorm(W,2,2)>1.5307740677166393);
g(tooLarge) = 1/vecnorm(W(tooLarge,:),2,2);
W = g.*W;
toReplace = find(sheepIsSeen==0);
if (~isempty(toReplace))
    meanDogsVel = mean(W(find(sheepIsSeen~=0),:),1);
    W(toReplace,1) = meanDogsVel(1);
    W(toReplace,2) = meanDogsVel(2);
end
U(1:Ndogs,:) = W;


%------------ Uncomment this and line 251 for stalking ------%
% indices = find(stalk);
% 
% stalkDir = CMs(indices,:)-X(indices,:);
% stalkDirNorm = vecnorm(stalkDir,2,2);
% if (stalkDirNorm~=0)
%     stalkDir = stalkDir./stalkDirNorm;
% end
% U(indices,:) = 0.02.*(stalkDir);
%------------------------------------------------------------%

output = {U, goalLocationArr,P,LastSeenCell,XmemCell,exponentialDecay,stopSimul,SubflockCells,areaOfAlpha/areaConvHull};