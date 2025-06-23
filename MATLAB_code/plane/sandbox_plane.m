%% Parameters
clear all;
N = 50;
Ndogs = 1;
L = 3.3;
nu = ones(N,1);
tmax = 1000; % set to 50 for shorter runtime
pressures = zeros(tmax,1);
CMDrift = zeros(tmax,1);

vMaxSheep = 1/2;
HomingDistance = 10*L;
%---------------------------
% syms x y
% scalarField(x,y) = 0;%5*exp((-x^2-y^2)/400);
% %scalarField(x,y) = x+y^2;
% v = [x,y];
% vectorField = gradient(scalarField,v);
% vecFunc = @(x,y) vectorField(x,y);
%---------------------------

filename='data.mat';
erase='';

% Display
display=true;   % if true, will plot agents when code is run
fixframe = true;
frameinrad = 30;

% no idea
fdim = 1;
M0 = L;
if fdim==2
    M0 = pi*L^2/2;
end



%% Initial conditions

ic_rad=sqrt(N)*L;

% sets seed
rng(4);
FourierCoeff = rand([10 10 3]);

% assign to X the value of ic_rad
X=zeros(N,2);
% multiplied by a random number drawn from the interval [-1,1]
X(1:Ndogs,:) = 4*ic_rad*(2*rand(Ndogs,2) - 1)-40;
X(Ndogs+1:N,:) = ic_rad*(2*rand(N-Ndogs,2) - 1)-10;

Nsheep = N-Ndogs;
rainbowPalette = hsv; % better for colourblind to use cool
colours = (mod(1:Nsheep, Nsheep)+1)';
angles = atan2(X(Ndogs+1:N,2),X(Ndogs+1:N,1));
[~,permutation] = sort(angles);
X(Ndogs+1:N,:) = circshift(X(permutation+Ndogs,:),(Nsheep-mod(Nsheep,2))/2);

% sets seed
rng(500);

sheepThetaVision = ((pi/180)*(306-191)).*rand(N,1) + (pi/180)*191;

ang = 2*pi*rand(N,1); % N x 1
U = [cos(ang) sin(ang)]; % N x 2, random starting orientations
U(1:Ndogs,:) = 0;

X_0 = X;
%% Initialize targets

% assigns subpopulations 1 and 2
% a N x 1 vector of all ones
ids = ones(N,1);
ids(1:N/2) = 2;



% this makes 3 targets:
% tar_rad = 4;
% tar_ang = 2*pi/3;
% tar = Target({tar_rad*[1 0],tar_rad*[cos(tar_ang) sin(tar_ang)],tar_rad*[cos(2*tar_ang) sin(2*tar_ang)]});
tar = Target();

% Our targets
%tar = Target([0 0]);
% tar = Target({[0 0],[20 20]});



dogTar = Target([20 12.5]);
%dogTar = Target();

%% Preallocation of some intermediate variables

% initialize all of these as zero vectors
[r,a,h0,h,s,l,rho] = deal(zeros(N,2));
prevSheepTars = zeros(Ndogs,1);

% initialize as N x 2 of zeroes
U1 = zeros(N,2);
Q = zeros(N,2);



%% Preallocation of time series variables
U_t = zeros(N,2,tmax);
DT_t = cell(tmax,1);
X_T = zeros(N,2,tmax);
Polarization_t = zeros(tmax,1);
Infections_t = cell(N-Ndogs,tmax);
spreadDivider = 1;
expDecay = zeros(N,2);
%enforce sheep speed limit
g = ones(N,1);
sizeOfVel = vecnorm(U,2,2);
tooLarge = find(sizeOfVel>vMaxSheep);
g(tooLarge) = vMaxSheep./sizeOfVel(tooLarge);
U(Ndogs+1:end,:) = g(Ndogs+1:end).*U(Ndogs+1:end,:);
LastSeen = zeros(N,1);
equil = 0;
TimeStratifiedX = 0;

scalarF = "fence";

%% Time for loop
for t = 1:tmax
    
    %---------------------------------%
 
    DT = delaunayTriangulation(X);
    [ConvexHull,hullArea] = convexHull(DT);
    [nbhd, nearest, d] = neighborhoods(DT);

    %---------------------------------%
    %keeping track of info at each time step
    DT_t{t} = DT;
    U_t(:,:,t) = U;
    X_T(:,:,t) = X;
    Polarization_t(t) = polarization(U,Ndogs);
    CMDrift(t) = vecnorm(mean(X(Ndogs+1:N,:)-[20 20]),2,2);
    pressures(t) = voronoiPressure(DT);

    %---------------------------------%

    % get alignment vector scaled by current velocity
    a = alignTo(X,U,nbhd,'expReciprocal', sheepThetaVision, Ndogs,vMaxSheep);

    %---------------------------------%
    prefVel = gradPreferenceField(X,FourierCoeff,scalarF);
    
    %Memory Dognamics
    if Ndogs>0

        dSV = dogSweepVoronoiMem2(X_T,U, DT, Ndogs, L, dogTar,t,LastSeen,scalarF,prefVel(1:Ndogs,:));
        U1 = dSV{1};
        equil = dSV{2};
        convHullIndices = dSV{3};
        LastSeen = dSV{4};
        TimeStratifiedX = dSV{5};
        plotTimeStratifiedX = true;
        expDecay = dSV{6};
        if (dSV{7})
            break
        end
    else
        convHullIndices = convexHull(DT);
        plotTimeStratifiedX = false;
    end
    %---------------------------------%
    %get the repulsion vector and value of sigma curve for each agent

    C = dogRepulsion(X, U1,DT, L, Ndogs,'expReciprocal','dogExpReciprocal',nbhd);
    r = C{1};
    s = C{2};
    nu = C{3}./2;

    %---------------------------------%
    %Get homing vector
    h0=homeToTarget(tar,X);

    %Turn off homing for agents too far away to feel the attraction
    hNorm = vecnorm(h0,2,2);
    indicesToZero = find(hNorm > HomingDistance);
    h0(indicesToZero,:) = 0;

    %normalize the homing vector with sigma
    h = (1-abs(s)) .* h0./vecnorm(h0,2,2);   % rescale to length 1-s
    h(isnan(h)) = 0;                    % protect 0 entries
    %---------------------------------%

    %add up all contributions to the velocity; divide by 5 for sheep:dog
    %speed ratios

    U1(Ndogs+1:N,:) = (r(Ndogs+1:N,:) + h(Ndogs+1:N,:) + nu(Ndogs+1:N).*a(Ndogs+1:N,:)+prefVel(Ndogs+1:N,:))./((1 + nu(Ndogs+1:N)));

    %---------------------------------%
    % calculate rho based on model choice
    [~,l] = voronoiProjectToBoundary(DT,U1);
    M = l;
    if fdim==2
        M = voronoiForwardArea(DT,U1);
    end
    %---------------------------------%

    % plot
    if display==true
        showAgents(X,U,tar,DT,fixframe,Ndogs,frameinrad,t,equil,convHullIndices,TimeStratifiedX,plotTimeStratifiedX,expDecay,colours,rainbowPalette,scalarF);
    end
    %---------------------------------%

    % Time step
    U = U1;
    U(Ndogs+1:N,:) = tanh(M(Ndogs+1:N,:)/M0).*U(Ndogs+1:N,:);

    %---------------------------------%


    %enforce cosmic sheep speed limit
    g = ones(N,1);
    sizeOfVel = vecnorm(U,2,2);
    tooLarge = find(sizeOfVel>vMaxSheep);
    g(tooLarge) = vMaxSheep./sizeOfVel(tooLarge);
    U(Ndogs+1:end,:) = g(Ndogs+1:end).*U(Ndogs+1:end,:);
    %---------------------------------%

    X = X + U;




    % display progress
    barwidth = 60;
    frac = floor(t/tmax * barwidth);
    fracpct = t/tmax*100;
    sofar = repmat(sprintf('\x0023'),1,frac);
    yet = repmat('-',1,barwidth-frac);
    progbar = ['[' sofar yet ']' ' ' sprintf('%04.1f\t%s',fracpct,filename) newline];
    fprintf([erase progbar]);
    erase = repmat(sprintf('\b'),1,length(progbar));


end
% 
S = unifMixMetric(N,Ndogs,delaunayTriangulation(X(Ndogs+1:end,:)),X_0);
%disp(S)

S1 = DiffOfExponentials(N,Ndogs,DT_t{1},DT_t{t});
disp(S1)
times = uint32(1):uint32(t);

fig2 = figure(2);
scatter(times, pressures(times));
title("Pressure vs Time")

fig3 = figure(3);
scatter(times, CMDrift(times));
title("Magnitude of Center of Mass Drift From (20,20) vs Time")

fig4 = figure(4);
plot(times,Polarization_t(times))
title("Polarization vs Time")

save("data.mat",'DT_t','U_t');