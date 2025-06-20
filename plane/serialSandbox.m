function serialSandbox(filename,N,L,nu,tmax,position_seed,angle_seed,Ndogs)
arguments
    filename (1,:) char
    N (1,1) double {mustBeInteger,mustBePositive}
    L (1,1) double {mustBeNonnegative}
    nu (1,1) double {mustBePositive}
    tmax (1,1) double {mustBeInteger,mustBePositive}
    position_seed (1,1) double {mustBeInteger} = 0
    angle_seed (1,1) double {mustBeInteger} = 0
    Ndogs (1,1) double {mustBeInteger} = 0
end

%% Parameters
erase='';

fdim = 1;
M0 = L;
if fdim==2
    M0 = pi*L^2/2;
end

%% Initial conditions
ic_rad=.5*sqrt(N*pi/4/.91);
rng(position_seed);
X = ic_rad*(2*rand(N,2) - 1);
rng(angle_seed);

sheepThetaVision = ((pi/180)*(306-191)).*rand(N,1) + (pi/180)*191;
ang = 2*pi*rand(N,1);
U = [cos(ang) sin(ang)];

%% Initialize targets
tar = Target([0 0]);
dogTar = Target([20 20]);

%% Preallocation of time series variables
U_t = zeros(N,2,tmax);
DT_t = cell(tmax,1);
U1 = zeros(N,2);

prevSheepTars = zeros(Ndogs,1);

for t = 1:tmax
    %% Compute step
    % make graph based on current positions of agents at time t
    DT = delaunayTriangulation(X);
    %[ConvexHull,hullArea] = convexHull(DT);
    % save in guy
    DT_t{t} = DT;

    % save current displacement vectors in guy at time t
    U_t(:,:,t) = U;

    [nbhd, nearest, d] = neighborhoods(DT);


    if Ndogs>0
        dSV = dogSweepVoronoi(X,U, DT, Ndogs, L, dogTar, nbhd,prevSheepTars);
        U1 = dSV{1};
        prevSheepTars = dSV{2};
        equil = dSV{3};
    end

    C = dogRepulsion(X, U1, L, nbhd, nearest, d, Ndogs,'expReciprocal','dogExpReciprocal');
    r = C{1};
    s = C{2};

    % alignment
    a = alignTo(X,U,nbhd,'expReciprocal', sheepThetaVision, Ndogs);

    % homing
    h0=homeToTarget(tar,X);
    h = (1-abs(s)) .* h0./vecnorm(h0,2,2);   % rescale to length 1-s
    h(isnan(h)) = 0;                    % protect 0 entries

    % direction
    U1(Ndogs+1:N,:) = (r(Ndogs+1:N,:) + h(Ndogs+1:N,:) + nu.*a(Ndogs+1:N,:))./(5*(1 + nu));

    % speed
    [~,l] = voronoiProjectToBoundary(DT,U1);
    M = l;
    if fdim==2
        M = voronoiForwardArea(DT,U1);
    end

    %% Update
    U = U1;
    U(Ndogs+1:N,:) = tanh(M(Ndogs+1:N,:)/M0).*U(Ndogs+1:N,:);
    X = X + U;

    %% display progress
    barwidth = 60;
    frac = floor(t/tmax * barwidth);
    fracpct = t/tmax*100;
    sofar = repmat(sprintf('\x0023'),1,frac);
    yet = repmat('-',1,barwidth-frac);
    progbar = ['[' sofar yet ']' ' ' sprintf('%04.1f\t%s',fracpct,filename) newline];
    fprintf([erase progbar]);
    erase = repmat(sprintf('\b'),1,length(progbar));


end

save(filename,'N','L','nu','Ndogs','position_seed','angle_seed','DT_t','U_t');
end