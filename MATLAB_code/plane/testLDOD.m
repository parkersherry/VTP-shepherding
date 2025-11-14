X = 1000*(2*rand(100,2) - 1);
rArr = linspace(1,10000,100);

dod = zeros(size(rArr));

for i=1:numel(rArr)
    dod(i) = LDOD(X,rArr(i));
end

scatter(rArr,dod)