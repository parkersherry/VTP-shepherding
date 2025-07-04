N=ceil(linspace(200,1000,20));

meanD = zeros(size(N));
minD = zeros(size(N));
maxD = zeros(size(N));

for k=1:numel(N)
    disp(k)
    temp = zeros(1,100);
    for i=1:100
        X = sqrt(N(k))*(2*rand(N(k),2) - 1);
        del = delaunay(X);
        g = graph(del, del(:, [2 3 1]));
        A = adjacency(g);
        A = A | A';
        g = graph(A);
        D = distances(g);
        %v = max(D,[],'all');
        temp(i) = max(D,[],'all');
    end
    meanD(k) = mean(temp);
    maxD(k) = max(temp);
    minD(k) = min(temp);
end
fig = figure();
hold on
scatter(N,meanD, 'filled');
% scatter(N,maxD, 'filled');
% scatter(N,minD, 'filled');
%fitobject = fit(N,meanD,)
%plot(N,log(N)./log(1.5))
disp([N;maxD])
xlabel('N');
ylabel('min max of D matrix')