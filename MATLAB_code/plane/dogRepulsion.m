function C = dogRepulsion(X,U,DT,L,Ndogs,transition_func, dog_trans_func,Firstnbhd)
[nbhd,nearest,d] = neighborhoods(DT,2);

fun = @(x) transition(x,transition_func);
funDog = @(x) transition(x,dog_trans_func);

%VERY NOT CONFIDENT IN THIS
dogL = 65*L/2;

LArrRepulsion = L.*ones(numel(d),1);

LArrAttractionToCM = L.*ones(numel(d),1);
ToCMStrength = ones(numel(d),1);
% logical array which returns 0 if dog is nearest neighbour
nearestIsDog = nearest >= Ndogs;

% agent-agent repulsion (scaled later)
r = X - X(nearest,:);   % vector from agent to nearest neighbour
rnorm = vecnorm(r,2,2); 
r = r./rnorm;
N = numel(X)/2;
% big dog repulsion
nu = ones(N,1);

rToDog = zeros(size(X));

for j = Ndogs+1:N

    allNbhd = [nbhd{j,1} nbhd{j,2}];
    allNbhd = unique(allNbhd);
    numJ = numel(allNbhd);

    for i = 1:numJ
        if allNbhd(i) <= Ndogs
            
            %repulsion strength is minimized when you're perp to the
            %direction of the dog's motion
            rToDog(j,:) = X(j,:) - X(allNbhd(i),:);
            rToDogNorm = vecnorm(rToDog(j,:),2,2);
            cosTheta = dot(rToDog(j,:)./rToDogNorm, U(allNbhd(i),:)./vecnorm(U(allNbhd(i),:),2,2));
            
            %sheep should not chase dogs
            if cosTheta <= 0
                cosTheta = 0;
            end

            sDog = funDog(rToDogNorm/dogL);
            rToDog(j,:) = 0.5 * (1+cosTheta)*sDog .* rToDog(j,:)./rToDogNorm;

            %values for nu, distance for repulsion and distance of CM
            %attraction for sheep near a dog
            nu(j) = 1/10;
            LArrRepulsion(j) = LArrRepulsion(j)/5;
            LArrAttractionToCM(j) = LArrAttractionToCM(j)*5;
            ToCMStrength(j) = 5;
        end
    end
end

s = arrayfun(fun, d./LArrRepulsion);
r = s .* r;  % rescale to length s
toCM = ToCMStrength.*goToCM(X,Firstnbhd,N,Ndogs,'linear2',LArrAttractionToCM);

% ignore flock repulsion if dog is nearest neigbhour
r(:,1) = r(:,1) .* nearestIsDog;
r(:,2) = r(:,2) .* nearestIsDog;

r = r + rToDog + toCM;
C = {r,s,nu};
end