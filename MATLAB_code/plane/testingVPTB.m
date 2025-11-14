
t = 1;
dt = 1/100;

X = [[1 1];[-1 1];[1 -1];[-1 -1];[0.5 0.5]];
U = [[1 0];[-1 0];[2 0];[0.5 -1];[-0.5 0.5]];
DT = delaunayTriangulation(X);
N = size(X(:,1));
N = N(1);
Ndogs=1;
disp(N)
tmax = 2;
X_T = zeros(N,2,2);
X_T(:,:,1) = X;
X_T(:,:,2) = X;

memDuration = 1;

tar = Target([]);
[nbhd,~,~] = neighborhoods(DT,1);
L = 3/4;
Larr = ones(N,1)*L;

sheepThetaVision = zeros(5,1);
Xmem = cell(Ndogs,1);
[Xmem{1:Ndogs}] = deal(zeros(N,2));

LastSeen = cell(Ndogs,1);
[LastSeen{1:Ndogs}] = deal(zeros(N,1));
out = dogMovementScheme(X_T, U, DT, Ndogs, L, tar, t, dt, LastSeen, Xmem,'zero','zero',Inf, memDuration);

disp(out{1})