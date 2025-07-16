function a = serialSandbox(filename,nuIn,N,L,alpha,tmax,dogTar,memDuration,position_seed,angle_seed,Ndogs,boundary,vMaxSheep)
arguments
    filename (1,:) char
    nuIn
    N (1,1) double {mustBeInteger,mustBePositive}
    L (1,1) double {mustBeNonnegative}
    alpha
    tmax (1,1) double {mustBeInteger,mustBePositive}
    dogTar
    memDuration
    position_seed (1,1) double {mustBeInteger} = 0
    angle_seed (1,1) double {mustBeInteger} = 0
    Ndogs (1,1) double {mustBeInteger} = 0
    boundary = 'zero'
    vMaxSheep = 1/2;
end

%% Parameters

Xmem = zeros(N,2);
HomingDistance = 10*L;

erase='';


% no idea
fdim = 1;
M0 = L;
if fdim==2
    M0 = pi*L^2/2;
end



%% Initial conditions
nu = ones(N,1);
Nsheep = N-Ndogs;

ic_radSheep=sqrt(Nsheep)*L;
ic_radDog=sqrt(Ndogs)*L;
% sets seed
rng(position_seed);
FourierCoeff = zeros([10 10 3]);

% assign to X the value of ic_rad
X=zeros(N,2);
% multiplied by a random number drawn from the interval [-1,1]
X(1:Ndogs,:) = ic_radDog*(2*rand(Ndogs,2) - 1)-sqrt(2*Nsheep)*L - sqrt(2*Ndogs)*L;
X(Ndogs+1:N,:) = ic_radSheep*(2*rand(Nsheep,2) - 1)-sqrt(Nsheep)*L;


% sets seed
rng(angle_seed);

sheepThetaVision = ((pi/180)*(306-191)).*rand(N,1) + (pi/180)*191;

ang = 2*pi*rand(N,1); % N x 1
U = [cos(ang) sin(ang)]; % N x 2, random starting orientations
U(1:Ndogs,:) = 0;

X_0 = X;
%% Initialize targets

tar = Target();


%% Preallocation of some intermediate variables


% initialize as N x 2 of zeroes
U1 = zeros(N,2);


%% Preallocation of time series variables
U_t = zeros(N,2,tmax);
DT_t = cell(tmax,1);
Polarization_t = zeros(1,tmax);
X_T = zeros(N,2,tmax);
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
DistMatrixCell = cell(1,tmax);
alphaHull_T = cell(1,tmax);
scalarF = boundary;

%% Time for loop
for t = 1:tmax
    
    %---------------------------------%
 
    DT = delaunayTriangulation(X);
    [nbhd, ~, ~] = neighborhoods(DT);

    %---------------------------------%
    %keeping track of info at each time step
    DT_t{t} = DT;
    U_t(:,:,t) = U;
    X_T(:,:,t) = X;
    % DistMatrixCell{t} = distanceMixMetric(N,Ndogs,X,X);
    Polarization_t(t) = polarization(U,Ndogs);
    %---------------------------------%

    % get alignment vector scaled by current velocity
    a = alignTo(X,U,nbhd,'expReciprocal', sheepThetaVision, Ndogs,vMaxSheep);

    %---------------------------------%
    prefVel = gradPreferenceField(X,FourierCoeff,scalarF);

    %Memory Dognamics
    if Ndogs>0

        DMS = dogMovementScheme(X_T,U, DT, Ndogs, L, dogTar,t,LastSeen,Xmem,scalarF,prefVel(1:Ndogs,:),1,memDuration);
        U1 = DMS{1};
        LastSeen = DMS{4};
        alphaHull_T{t} = DMS{3};
        Xmem = DMS{5};

        if (DMS{7})
            break
        end
    end
    %---------------------------------%
    %get the repulsion vector and value of sigma curve for each agent

    C = sheepMovementScheme(X, U1,DT, L, Ndogs,'expReciprocal','dogExpReciprocal',nbhd,tar,HomingDistance);
    r = C{1};
    h = C{2};
    nu = C{3}.*nuIn;


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


    % Time step
    U1(Ndogs+1:N,:) = tanh(M(Ndogs+1:N,:)/M0).*U1(Ndogs+1:N,:);

    %---------------------------------%


    %enforce cosmic sheep speed limit
    g = ones(N,1);
    sizeOfVel = vecnorm(U1,2,2);
    tooLarge = find(sizeOfVel>vMaxSheep);
    g(tooLarge) = vMaxSheep./sizeOfVel(tooLarge);
    U(Ndogs+1:end,:) = g(Ndogs+1:end).*U1(Ndogs+1:end,:);
    %---------------------------------%

    U(1:Ndogs,:) = U1(1:Ndogs,:);

    X = X + U;
    %showAgents(X,U,tar,DT,false,Ndogs,20,t,[0 0],alphaHull,TimeStratifiedX,false,expDecay,0,hsv,'zero');

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

% times = linspace(1,tmax,tmax);
% modfun = @(a,t) a(1).*tanh(a(2).*t./log(2));
% mdl = fitnlm(times,muMixing,modfun,[6 3]);
% a = mdl.Coefficients.Estimate;
save(filename,'tmax','position_seed','angle_seed','Polarization_t');
%save(filename,'position_seed','angle_seed','N','aEstimated');