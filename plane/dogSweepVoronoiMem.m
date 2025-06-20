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

function output = dogSweepVoronoiMem(X_T, U,DT, Ndogs, L, tar,previousSheepTars ,t ,LastSeen,spreadDivider)

arguments
    X_T
    U
    DT
    Ndogs % int8 {mustBePositive}
    L double {mustBePositive}
    tar
    previousSheepTars
    t
    LastSeen
    spreadDivider
end
% toggles whether to prioritize being collinear with
% local CM and dog target
Shepherd = true;

X = X_T(:,:,t);

N = numel(X)/2;
Xsheep = X(Ndogs+1:end,:);

recallDelay = 240;

[nbhd, nearest, ~] = neighborhoods(DT,2);
DT_sheep = delaunayTriangulation(X(Ndogs+1:N,:));
[C, ~] = convexHull(DT_sheep);
chullTar = Target(Xsheep(C,:));

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

    % calculate local CM
    J = numel(SheepNbhdPast);
    for j=1:J
        X_j = X_T(:,:,LastSeen(SheepNbhdPast(j)));
        TotalX(SheepNbhdPast(j),:) = X_j(SheepNbhdPast(j),:);
    end

    %-----------------------%
    %if you can see the residual image of the sheep then obvi there is no
    %sheep there and you should forget it
    ToDel = zeros(size(SheepNbhdPast));
    inside = inhull(TotalX(SheepNbhdPast,:),X(C,:));
    for j=1:numel(SheepNbhdPast)
        if (inside(j)==1)
            continue
        end
        Xcoords = [X(i,1) TotalX(SheepNbhdPast(j),1)];
        %disp(size(Xcoords))
        Ycoords = [X(i,2) TotalX(SheepNbhdPast(j),2)];
        [xIntersect,~] = polyxpoly(Xcoords,Ycoords,X(C,1),X(C,2));
        if (isempty(xIntersect))
            ToDel(j) = 1;
        end
    end
    ToDel = find(ToDel);
    TotalX(SheepNbhdPast(ToDel),:) = 0;
    LastSeen(SheepNbhdPast(ToDel)) = 0;

    SheepNbhdPast(ToDel) = [];
    %-----------------------%
    CurrAllNbhd = [nbhd{i,1} nbhd{i,2}];

    CurrAllNbhd = unique(CurrAllNbhd);
    CurrSheepnbhd = CurrAllNbhd(CurrAllNbhd>Ndogs);
    LastSeen(CurrSheepnbhd) = t;
    TotalX(CurrSheepnbhd,:) = X(CurrSheepnbhd,:);
    SheepNbhd = unique(union(SheepNbhdPast,CurrSheepnbhd));

    exponentialDecay = ExpAmplitude.*exp(decayRate.*(LastSeen-t)./recallDelay) + ExpShift;
    CM = mean(TotalX(SheepNbhd,:).*exponentialDecay(SheepNbhd));

    equilibrium = findDogEquil(CM,L , tar,numel(SheepNbhd));

    % finds sheep agent with indexMax

    xDog = X(i,:);
    CMtoDog = xDog - CM;
    CMtoDog = CMtoDog ./ vecnorm(CMtoDog,2,2);

    DogToSheep = xDog-TotalX(SheepNbhd,:);
    DogToSheep = DogToSheep./vecnorm(DogToSheep,2,2);

    thetaArray = acos(DogToSheep*transpose(CMtoDog));
    % toBias = find(SheepNbhd == previousSheepTars(i));
    % bias = 0;
    % if ~isempty(toBias)
    %     thetaArray(toBias) = thetaArray(toBias) + bias;
    % end
    [~,~,intersectionIndices] = intersect(SheepToCollectArr,SheepNbhd);
    thetaArray(SheepNbhd(intersectionIndices)) = 0;
    [~,tempIndex] = max(thetaArray);
    SheepToCollect = SheepNbhd(tempIndex);
    SheepToCollectArr(i) = SheepToCollect;


    h = homeToTarget(chullTar,X(i,:),1,false);
    distToConvHull = vecnorm(h,2,2);
    if (distToConvHull ~= 0)
        h = h./distToConvHull;
    end

    if Shepherd
        % equilibrium ON protocol
        perp = transpose([CMtoDog(2) (-1*CMtoDog(1))]);
        Xprojected = TotalX(SheepNbhd,:)*perp;
        maxSpread = max(Xprojected);
        minSpread = min(Xprojected);
        spread = abs(maxSpread-minSpread);

        if spread > sqrt(numel(SheepNbhd))*L/spreadDivider
            

            spreadDivider = 5/4;

            XsheepMaxIndex = TotalX(SheepToCollect,:);

            goalLocation = (XsheepMaxIndex - CM);
            goalLocation =  goalLocation./vecnorm(goalLocation,2,2);
            goalLocation = goalLocation.*L;
            goalLocation = goalLocation + XsheepMaxIndex;
            dogToGoal = goalLocation - xDog;
            
            if vecnorm(dogToGoal,2,2) < dogStalkLength
                stalk(i) = 1;
                continue
            end

            dogToGoal = dogToGoal ./ vecnorm(dogToGoal,2,2);
            dogS = 0;

            if nearest(i)<Ndogs+1
                dogS = transition(vecnorm(xDog - X(nearest(i),:),2,2)/(L), 'expReciprocal');
            end

            s = transition(distToConvHull/L, 'expReciprocal');
            dog_dogRepel = (xDog - X(nearest(i),:))./vecnorm(xDog - X(nearest(i),:),2,2);
            DogDisplacementVec(i,:) = (-1*repelConvHullStrength*s).* h + (1-s) .* dogToGoal + (dogS).*dog_dogRepel;
            goalLocationArr(i,:) = goalLocation;
        else
            spreadDivider = 1;

            % vector from dog to equilibrium
            DogToEquil = equilibrium - xDog;
            DogToEquilNorm = vecnorm(DogToEquil,2,2);
            if DogToEquilNorm < dogStalkLength
                stalk(i) = 1;
                continue
            end
            sEquil = transition(vecnorm((xDog - equilibrium)/(L),2,2), 'expReciprocal');
            s = transition(distToConvHull/L, 'expReciprocal');

            DogDisplacementVec(i,:) = (-1*repelConvHullStrength*s).* h + (1-abs(sEquil)) .* DogToEquil ./ DogToEquilNorm;
            goalLocationArr(i,:) = equilibrium;

        end
    else
        % see above sans equilibrium
        XsheepMaxIndex = TotalX(SheepToCollect,:);

        goalLocation = (XsheepMaxIndex - CM);
        goalLocation =  goalLocation./vecnorm(goalLocation,2,2);
        goalLocation = goalLocation.*L;
        goalLocation = goalLocation + XsheepMaxIndex;
        dogToGoal = goalLocation - xDog;
        if vecnorm(dogToGoal,2,2) < dogStalkLength
            stalk(i) = 1;
            continue
        end
        dogToGoal = dogToGoal ./ vecnorm(dogToGoal,2,2);
        dogS = 0;
        if nearest(i)<Ndogs+1
            dogS = transition(vecnorm(xDog - X(nearest(i),:),2,2)/(L), 'expReciprocal');
        end
        s = transition(distToConvHull/L, 'expReciprocal');
        dog_dogRepel = (xDog - X(nearest(i),:))./vecnorm(xDog - X(nearest(i),:),2,2);
        DogDisplacementVec(i,:) = (-1*repelConvHullStrength*s).* h + (1-s) .* dogToGoal + (dogS).*dog_dogRepel;
        goalLocationArr(i,:) = goalLocation;
    end
    CMs(i,:) = CM;
    %Xtotals{i} = TotalX;

end

%this 3 is arbitrary
clear("chullTar")
W = U(1:Ndogs,:)+(DogDisplacementVec);
g = ones(Ndogs,1);
tooLarge = find(vecnorm(W,2,2)>0.75);
g(tooLarge) = 1/vecnorm(W(tooLarge,:),2,2);

U(1:Ndogs,:) = g.*W;
indices = find(stalk);

stalkDir = CMs(indices,:)-X(indices,:);
stalkDir = stalkDir./vecnorm(stalkDir,2,2);
U(indices,:) = 0.02.*(stalkDir);

output = {U,SheepToCollectArr, goalLocationArr,C,LastSeen,TotalX,spreadDivider,exponentialDecay};