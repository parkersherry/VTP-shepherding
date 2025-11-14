function DOD = LDOD(X,r)

DOD = 0;
N = size(X,1);

X(end+1,:) = [10^6,10^6];
X(end+1,:) = [-10^6,10^6];
X(end+1,:) = [-10^6,-10^6];
X(end+1,:) = [10^6,-10^6];
DT = delaunayTriangulation(X);
[V,C] = voronoiDiagram(DT);
disp(C)
n=100;
theta = linspace(0,2*pi,n);
for i = 1:N
    cell = C{i};
    verts = V(cell,:);
    vRegion = polyshape(verts);

    cx = X(i,1);
    cy = X(i,2);
    x = cx + r*cos(theta);
    y = cy + r*sin(theta);
    circ = polyshape(x,y);

    limited = intersect(circ,vRegion);
    [v1,v2] = boundary(limited);

    A = poly_area(v1,v2);
    if isnan(A)
        A = Inf;
    end
    DOD = DOD + A;
    
    
    
end

DOD = DOD/(N*pi*r*r);
