warning('off','all')
theta = linspace(0,2*pi,50);
r = 1;
X =zeros(100,2);
for i=1:50
    [X(i,1), X(i,2)] = pol2cart(theta(i),r);
    [X(i+50,1), X(i+50,2)] = pol2cart(theta(i),r);
end
X(51:end,:) = X(51:end,:)+20;

subflocksCell = {};

DT = delaunayTriangulation(X);
% scatter(X(:,1),X(:,2),'filled')
hold on
% triplot(DT)
S = getSubflocks(X, X, 1, 100, 0, subflocksCell, 0);
disp(S)
for i=1:numel(S)
    scatter(S{i}(:,1), S{i}(:,2), 'filled')
end
