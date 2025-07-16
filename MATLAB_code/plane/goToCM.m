function G = goToCM(X, Firstnbhd, N, Ndogs,trans_func,LArr,L,cutoff)

arguments

    X   (:,2) double
    Firstnbhd 
    N   double  {mustBePositive}
    Ndogs   double
    trans_func 
    LArr
    L
    cutoff = 1000*L
end
G = zeros(N,2);
for i=1:N
    if i<Ndogs+1
        continue
    end
    FirstSheepNbhd = Firstnbhd{i}(Firstnbhd{i}>Ndogs);
    Xnbhd = X(FirstSheepNbhd,:);
    distForCutoff = vecnorm(Xnbhd-X(i),2,2);
    XCutoff = Xnbhd(find(distForCutoff<cutoff),:);
    if (~isempty(XCutoff))
        CM = mean(XCutoff,1);
        G(i,:) = CM - X(i,:);
    end
end
Gnormed = vecnorm(G,2,2);
f = @(x) transition(x,trans_func);
s = arrayfun(f,Gnormed./LArr);

G = s.*G;

end