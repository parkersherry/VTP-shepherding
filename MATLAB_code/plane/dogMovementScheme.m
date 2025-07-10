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
function output = dogMovementScheme(X_T, U,DT, Ndogs, L, tar ,t ,LastSeen,Xmem,scalarField,prefVel,alpha,memDuration)

arguments
    X_T
    U
    DT
    Ndogs % int8 {mustBePositive}
    L double {mustBePositive}
    tar
    t
    LastSeen
    Xmem
    scalarField
    prefVel
    alpha
    memDuration = 240
end
% toggles whether to prioritize being collinear with
% local CM and dog target
stopSimul = false;
X = X_T(:,:,t);

N = numel(X)/2;
Xsheep = X(Ndogs+1:end,:);

recallDelay = memDuration;

[nbhd, nearest, ~] = neighborhoods(DT,2);

shp = alphaShape(Xsheep);
alpha = criticalAlpha(shp,'one-region');
shp.Alpha = alpha;
[~,P] = boundaryFacets(shp);

if (~isempty(P))
    P(end+1,:) = P(1,:);
end
indicesConvHull = convhull(Xsheep);
ConvHullTarget = Target(Xsheep(indicesConvHull,:));
CMs = zeros(Ndogs,2);

% keeps track of sheep agent indices while finding maximizing agent
SheepToCollectArr = zeros(Ndogs,1);

% initialize dog agent displacement vector
DogDisplacementVec = zeros(Ndogs, 2);
Xtotals = cell(Ndogs);

% initialize positions to which dog agents move
goalLocationArr = DogDisplacementVec;

stalk = zeros(Ndogs,1);
dogStalkLength = 3.5;
%%%%%%%%%%THIS ONLY WORKS FOR ONE DOG
TotalX = zeros(N,2);

exponentialDecay = zeros(N,2);
decayRate = recallDelay*0.02;
ExpAmplitude = 0.43;
ExpShift = 0.48;

repelConvHullStrength = 1;
for i = 1:Ndogs


    if t < recallDelay
        timeThreshold = 0;
    else
        timeThreshold = t - recallDelay;
    end

    SheepNbhdPast = find(LastSeen > timeThreshold);
    forgotten = find((~(LastSeen > timeThreshold)).*(LastSeen~=0));
    

    % for j=1:J
    %     X_j = X_T(:,:,LastSeen(SheepNbhdPast(j)));
    %     TotalX(SheepNbhdPast(j),:) = X_j(SheepNbhdPast(j),:);
    % end

    %-----------------------%
    %if you can see the residual image of the sheep then obvi there is no
    %sheep there and you should forget it
    if (~isempty(SheepNbhdPast))
        ToDel = zeros(size(SheepNbhdPast));
        inside = inhull(Xmem(SheepNbhdPast,:),P);
        for j=1:numel(SheepNbhdPast)
            if (inside(j)==1)
                continue
            end
            Xcoords = [X(i,1) Xmem(SheepNbhdPast(j),1)];
            %disp(size(Xcoords))
            Ycoords = [X(i,2) Xmem(SheepNbhdPast(j),2)];
            [xIntersect,~] = polyxpoly(Xcoords,Ycoords,P(:,1),P(:,2));
            if (isempty(xIntersect))
                ToDel(j) = 1;
            end
        end
        ToDel = find(ToDel);
        Xmem(SheepNbhdPast(ToDel),:) = 0;
        LastSeen(SheepNbhdPast(ToDel)) = 0;
    
        SheepNbhdPast(ToDel) = [];
    end
    
    %-----------------------%
    CurrAllNbhd = [nbhd{i,1} nbhd{i,2}];

    CurrAllNbhd = unique(CurrAllNbhd);
    forgotten = setdiff(forgotten,CurrAllNbhd);
    CurrSheepnbhd = CurrAllNbhd(find(CurrAllNbhd>Ndogs));

    Xmem(CurrSheepnbhd,:) = X(CurrSheepnbhd,:);
    
    LastSeen(CurrSheepnbhd) = t;
    currMeanVel = mean(U(CurrSheepnbhd,:),1);

    %TotalX(CurrSheepnbhd,:) = X(CurrSheepnbhd,:);
    SheepNbhd = unique(union(SheepNbhdPast,CurrSheepnbhd));
    
    Xmem(forgotten,:) = 0;
    if strcmp(scalarField,'fence')
        indices = find(~((Xmem(SheepNbhd,1)<20).*(Xmem(SheepNbhd,2)<20)));
        SheepNbhd(indices) = [];
        if (isempty(SheepNbhd))
            disp("all sheep are outside the fence")
            stopSimul = true;
            break
        end
    end
    
    SheepNbhd( find(~any(Xmem(SheepNbhd,:)')), : ) = [];

    if (isempty(SheepNbhd))
        disp("The dog sees no sheep")
        stopSimul = true;
        break
    end

    exponentialDecay = ExpAmplitude.*exp(decayRate.*(LastSeen(SheepNbhd)-t)./recallDelay) + ExpShift;
    sumExpDecay=(sum(exponentialDecay));
    CM = sum(Xmem(SheepNbhd,:).*exponentialDecay,1)./sumExpDecay;
    equilibrium = findDogEquil(CM,L , tar,numel(SheepNbhd));

    % finds sheep agent with indexMax

    xDog = X(i,:);
    CMtoDog = xDog - CM;
    CMtoDogNorm = vecnorm(CMtoDog,2,2);
    if (CMtoDogNorm~=0)
        CMtoDog = CMtoDog ./CMtoDogNorm ;
    end

    DogToSheep = xDog-Xmem(SheepNbhd,:);
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


    XsheepMaxIndex = Xmem(SheepToCollect,:);

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
    dogS = 0;

    if nearest(i)<Ndogs+1
        dogS = transition(vecnorm(xDog - X(nearest(i),:),2,2)/(L), 'expReciprocal');
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
    dog_dogRepel = (xDog - X(nearestIndex,:));

    dog_dogRepelNorm = vecnorm(dog_dogRepel,2,2);
    if (dog_dogRepelNorm~=0)
        dog_dogRepel = dog_dogRepel./dog_dogRepelNorm;
    end
    DogDisplacementVec(i,:) =  (1-s) .* dogToGoal.*(1-sDriving)+ DogToEquil.*(sDriving);
    DogDisplacementVec(i,:) = DogDisplacementVec(i,:) -(1*repelConvHullStrength*s).* h + (dogS).*dog_dogRepel;
    goalLocationArr(i,:) = goalLocation;

    CMs(i,:) = CM;
    %Xtotals{i} = Xmem;
    Xmem(SheepNbhd,:) = Xmem(SheepNbhd,:) + currMeanVel;
end

if (strcmp(scalarField,'fence'))
    U(1:Ndogs,:) = U(1:Ndogs,:)+prefVel;
end


clear("chullTar")
W = U(1:Ndogs,:)+(DogDisplacementVec);
g = ones(Ndogs,1);
tooLarge = find(vecnorm(W,2,2)>1);
g(tooLarge) = 1/vecnorm(W(tooLarge,:),2,2);
U(1:Ndogs,:) = g.*W;
disp(find(isnan(U)))

indices = find(stalk);

stalkDir = CMs(indices,:)-X(indices,:);
stalkDirNorm = vecnorm(stalkDir,2,2);
if (stalkDirNorm~=0)
    stalkDir = stalkDir./stalkDirNorm;
end
U(indices,:) = 0.02.*(stalkDir);
output = {U, goalLocationArr,P,LastSeen,Xmem,exponentialDecay,stopSimul};