function S = DiffOfExponentials(N, Ndogs,DT_0, DT_tmax)

X_0 = DT_0.Points;
X_0 = X_0(Ndogs+1:end,:);
X_tmax = DT_tmax.Points;
X_tmax = X_tmax(Ndogs+1:end,:);

DT_0 = delaunayTriangulation(X_0);
DT_tmax = delaunayTriangulation(X_tmax);

E0 = DT_0.edges;
Etmax = DT_tmax.edges;

G0 = graph(E0(:,1),E0(:,2));
Gtmax = graph(Etmax(:,1),Etmax(:,2));

A0 = adjacency(G0);
Atmax = adjacency(Gtmax);

diff = expm(A0)-expm(Atmax);
diff = diff .* diff';
S = sum(diff,"all")/(N^2);