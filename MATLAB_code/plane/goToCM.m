function G = goToCM(X, Firstnbhd, N, Ndogs,trans_func,L,cutoff)

arguments

    X   (:,2) double
    Firstnbhd 
    N   double  {mustBePositive}
    Ndogs   double
    trans_func 
    L
    cutoff = 1000*L
end

G = zeros(N,2);
for i=1:N
    if i<Ndogs+1
        continue
    end
    FirstSheepNbhd = Firstnbhd{i}(Firstnbhd{i}>Ndogs);
    CM = mean(X(FirstSheepNbhd,:));
    G(i,:) = CM - X(i,:);
end
Gnormed = vecnorm(G,2,2);
f = @(x) transition(x,trans_func);
s = arrayfun(f,Gnormed./L);
s(find(Gnormed>cutoff)) = 0;
G = s.*G;

end