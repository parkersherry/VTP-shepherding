clear all
N=ceil(linspace(100,1000,5));

meanDiff = zeros(size(N));
minDiff = zeros(size(N));
maxDiff = zeros(size(N));

for k=1:numel(N)
    disp(k)
    temp =zeros(1,1000);
    for i=1:1000
        X1 = sqrt(N(k))*(2*rand(N(k),2) - 1);
        X2 = sqrt(N(k))*(2*rand(N(k),2) - 1);
        %S = DiffOfExponentials(N(k),0,delaunayTriangulation(X1),delaunayTriangulation(X2));
        S = distanceMixMetric(N(k),0,(X1),X2);
        temp(i)=S;
    end
    maxDiff(k)=max(temp);
    minDiff(k)=min(temp);
    meanDiff(k)=mean(temp);
end
fig3 = figure(3);
hold on
scatter(N,maxDiff,'filled','black')
scatter(N,meanDiff,'filled','red')
scatter(N,minDiff,'filled','blue')

title("NumEdges vs N")