function a = alignTo(X,U,nbhd,transition_func, sheepThetaVision,Ndogs,vMaxSheep)

arguments
    X (:,2) double
    U (:,2) double
    nbhd (:,1) cell {mustBeEqualSize(nbhd,U)}
    transition_func (1,:) char
    sheepThetaVision
    Ndogs
    vMaxSheep
end
adjust = false;

N = size(U,1);
velBaseds = ones(N,1);
a = zeros(N,2);
g = @(x) transition(x,transition_func);
w = @(x) g((1/pi)*acos(x));
for i=Ndogs+1:N
    ui = U(i,:);
    ui = ui/norm(ui,2);
    ui(isnan(ui)) = 0;
    J = nbhd{i};
    velBaseds(i) = 2.-g(vecnorm(max(U(J,:)))/vMaxSheep);
    for j=1:size(J,2)
        s_x = X(i,:) - X(nbhd{i}(j),:);
        s_x = s_x./vecnorm(s_x,2,2);
        if acos(dot(s_x, ui)) > sheepThetaVision(i)/2
            continue
        end
        uj = U(J(j),:);
        uj = uj/norm(uj,2);
        uj(isnan(uj)) = 0;
        s = dot(ui,uj);
        % due to numerical error, it can happen that |s|>1, if this is not
        % corrected, w(s) (which calls acos) will attempt to input an
        % imaginary value to g. The complex-valued acos function is
        % continuous at ±1 so this correction is well-behaved.
        if abs(s) > 1
            s = sign(s);
        end
        a(i,:) = a(i,:) + w(s)*uj;
    end
end
a = (1/6)*a;

s = ones(N,1);

s(Ndogs+1:end,:) = 2-arrayfun(g,abs(vecnorm(U(Ndogs+1:end,:),2,2)./vMaxSheep));
if adjust
    a = s.*a;
    a = velBaseds.*a;
end
end

% validation
function mustBeEqualSize(a,b)
if size(a,1) ~= size(b,1)
    eid = 'Size:equalSize';
    msg = 'second input must have same length as first';
    throwAsCaller(MException(eid,msg));
end
end