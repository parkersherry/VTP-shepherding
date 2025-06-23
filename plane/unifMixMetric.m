
function S = unifMixMetric(N, Ndogs, DT_final, X_0)

arguments
    N double {mustBePositive}
    Ndogs %double {mustBePositive}
    DT_final
    X_0
end

Nsheep = N - Ndogs;
del = delaunay(X_0(Ndogs+1:N,:));

g = graph(del, del(:, [2 3 1]));
A = adjacency(g);
A = A | A';
g = graph(A);

d = distances(g);
K = max(d,[],"all");
s = zeros(N,1);
[nbhd, nearest, dist] = neighborhoods(DT_final,1);


for i = 1:Nsheep
    J = numel(nbhd{i});
    dist = zeros(J,1);
    for j = 1:J
        dist(j) = d(i, nbhd{i}(j));
    end
    dist = dist - ones(J,1);
    s(i) = sum(dist);
    s(i) = s(i)/J;
end

S = 1/N * 1/(K-1) * sum(s);

end