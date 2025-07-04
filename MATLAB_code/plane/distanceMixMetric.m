function S = distanceMixMetric(N, Ndogs, X_0, X_f)

arguments
    N double {mustBePositive}
    Ndogs %double {mustBePositive}
    X_0
    X_f
end

Nsheep = N - Ndogs;
del0 = delaunay(X_0(Ndogs+1:N,:));

g = graph(del0, del0(:, [2 3 1]));
A0 = adjacency(g);
A0 = A0 | A0';
g0 = graph(A0);
D0 = distances(g0);

% delf = delaunay(X_f(Ndogs+1:N,:));
% 
% g = graph(delf, delf(:, [2 3 1]));
% Af = adjacency(g);
% Af = Af | Af';
% gf = graph(Af);
% Df = distances(gf);

% D0 = D0./mean(D0,'all');
% Df = Df./mean(Df,'all');
%S = sum(abs(D0-Df),'all')/(N-1)^2;
%S=sum((D0-Df).^2,'all')/(Nsheep^2-Nsheep);
S = D0;
