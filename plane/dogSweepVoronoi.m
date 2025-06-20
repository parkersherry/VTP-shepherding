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

function output = dogSweepVoronoi(X, U,DT, Ndogs, L, tar, firstNeighbors,previousSheepTars)

arguments
    X
    U
    DT
    Ndogs % int8 {mustBePositive}
    L double {mustBePositive}
    tar
    firstNeighbors
    previousSheepTars

end
% toggles whether to prioritize being collinear with
% local CM and dog target
shouldEqual = true;

N = numel(X)/2;
Nsheep = N - Ndogs;
Xsheep = X(Ndogs+1:end,:);

%dog repulsion range
Ldog = 65*L/2;

[nbhd, nearest, dist] = neighborhoods(DT,2);
DT_sheep = delaunayTriangulation(X(Ndogs+1:N,:));
[C, sheepArea] = convexHull(DT_sheep);
chullTar = Target(Xsheep(C,:));
CMs = zeros(Ndogs,2);
% for each dog, there is a vector filled with the indices of agents who are
% either first neighbours of the dog, or second neighbours on the convHull
I = cell(Ndogs,1);

% keeps track of sheep agent indices while finding maximizing agent
indexMax = zeros(Ndogs,1);

% initialize dog agent displacement vector
P = zeros(Ndogs, 2);

% initialize positions to which dog agents move
positions = P;

stalk = zeros(N,1);

LocalCM = zeros(N,2);

repulsionStrength = 1;

for i = 1:Ndogs
    allNbhd = [nbhd{i,1} nbhd{i,2}];
    h = homeToTarget(chullTar,X(i,:),1,false);
    distToConvHull = vecnorm(h,2,2);
    if (distToConvHull ~= 0)
        h = h./distToConvHull;
    end
    %%%IF h= [0 0] THIS CRASHES!!!

    dognbhd = allNbhd(allNbhd<Ndogs+1);
    sheepnbhd = allNbhd(allNbhd>Ndogs);

    % calculate local CM
    CM = mean(X(sheepnbhd,:));
    LocalCM(i,:) = CM;


    % ensures that members of I{i} are in sheepnbhd and on convHull or
    % first neighbours
    J = numel(nbhd{i,2});
    I{i} = nbhd{i,1};
    for j = 1:J
        if (ismember(nbhd{i,2}(j), C) && nbhd{i,2}(j)>Ndogs)
            I{i}(end+1) = sheepnbhd(j);
        end
    end
    equilibrium = findDogEquil(CM,L , tar,numel(I{i}));

    % finds sheep agent with indexMax
    MaxIndex = i;
    thetaMax = -1;
    K = numel(I{i});
    xDog = X(i,:);
    CMtoDog = xDog - CM;
    CMtoDog = CMtoDog ./ vecnorm(CMtoDog,2,2);
    for k = 1:K
        if I{i}(k) < Ndogs+1
            continue
        end
        if ismember(I{i}(k),indexMax(dognbhd))
            continue
        end
        xK = X(I{i}(k),:);
        bias = 0;
        if (I{i}(k)) == previousSheepTars(i)
            bias = 0; %%%BIAS HERE
        end
        V = xDog - xK;
        V = V ./ vecnorm(V,2,2);
        thetaK = acos(dot(CMtoDog,V)) + bias;
        if thetaMax < thetaK
            thetaMax = thetaK;
            MaxIndex = I{i}(k);
        end
    end

    indexMax(i) = MaxIndex;

    a = 0.02;

    if shouldEqual
        % equilibrium ON protocol
        perp = transpose([CMtoDog(2) (-1*CMtoDog(1))]);
        Xprojected = X(allNbhd,:)*perp;
        maxSpread = max(Xprojected);
        minSpread = min(Xprojected);
        spread = abs(maxSpread-minSpread);
        %disp(spread/(sqrt(Nsheep)*L))
        if spread > sqrt(numel(allNbhd))*L

            XsheepMaxIndex = X(MaxIndex,:);

            pos = (XsheepMaxIndex - CM);
            pos =  pos./vecnorm(pos,2,2);
            pos = pos.*L;
            pos = pos + XsheepMaxIndex;
            dogToPos = pos - xDog;
            if vecnorm(dogToPos,2,2) < 3.5
                stalk(i) = 1;
                % P(i,:) = -a .* CMtoDog;
                continue
            end
            dogToPos = dogToPos ./ vecnorm(dogToPos,2,2);
            dogS = 0;
            repulsionCM = CMtoDog;
            if nearest(i)<Ndogs+1
                dogS = transition(vecnorm(xDog - X(nearest(i),:),2,2)/(L), 'expReciprocal');
            end
            s = transition(distToConvHull/L, 'expReciprocal');
            dog_dogRepel = (xDog - X(nearest(i),:))./vecnorm(xDog - X(nearest(i),:),2,2);
            P(i,:) = (-1*repulsionStrength*s).* h + (1-s) .* dogToPos + (dogS).*dog_dogRepel;
            positions(i,:) = pos;

        else
            % vector from dog to equilibrium
            chase = equilibrium - xDog;
            chaseNorm = vecnorm(chase,2,2);
            if chaseNorm < 3.5
                stalk(i) = 1;
                % P(i,:) = -a .* CMtoDog;
                continue
            end
            sEquil = transition(vecnorm((xDog - equilibrium)/(L),2,2), 'expReciprocal');
            s = transition(distToConvHull/L, 'expReciprocal');

            P(i,:) = (-1*repulsionStrength*s).* h + (1-abs(sEquil)) .* chase ./ chaseNorm;
            positions(i,:) = equilibrium;
        end
    else
        % see above sans equilibrium
        XsheepMaxIndex = X(MaxIndex,:);

        pos = (XsheepMaxIndex - CM);
        pos =  2 .* pos./vecnorm(pos,2,2);
        %pos = pos.*R;
        pos = pos + XsheepMaxIndex;
        dogToPos = pos - xDog;
        if vecnorm(dogToPos,2,2) < 3.5
            stalk(i) = 1;
            % P(i,:) = -a .* CMtoDog;
            continue
        end
        dogToPos = dogToPos ./ vecnorm(dogToPos,2,2);
        dogS = 0;
        repulsionCM = CMtoDog;
        if nearest(i)<Ndogs+1
            dogS = transition(vecnorm(xDog - X(nearest(i),:),2,2)/(L), 'expReciprocal');
        end
        s = transition(distToConvHull/L, 'expReciprocal');
        dog_dogRepel = (xDog - X(nearest(i),:))./vecnorm(xDog - X(nearest(i),:),2,2);
        P(i,:) = (-1*repulsionStrength*s).* h + (1-s) .* dogToPos + (dogS).*dog_dogRepel;
        positions(i,:) = pos;
    end
    CMs(i,:) = CM;
end

%this 3 is arbitrary
clear("chullTar")
W = U(1:Ndogs,:)+(P./3);
g = ones(Ndogs,1);
tooLarge = find(vecnorm(W,2,2)>0.75);
g(tooLarge) = 1/vecnorm(W(tooLarge,:),2,2);

U(1:Ndogs,:) = g.*W;
indices = find(stalk);

stalkDir = LocalCM(indices,:)-X(indices,:);
stalkDir = stalkDir./vecnorm(stalkDir,2,2);
U(indices,:) = 0.02.*(stalkDir);

output = {U,indexMax, positions,C};