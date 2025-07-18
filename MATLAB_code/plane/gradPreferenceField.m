function grad = gradPreferenceField(X,c,scalarField)
%-----------------%
%inputs:
% -> X = the positions at which to find the gradient
% -> c = an array of coefficients for the Fourier components (to be
%        randomly generated
% -> scalarField = a string corresponding to the scalar field we want the
%                  gradient of

%outputs:
% -> grad = an array of size X, indicating the direction each element of X
%           should go in to follow the grad of the scalar field
%-----------------%
grad = zeros(size(X));
switch scalarField
    case '2gaussiansPositive'
        offset = 20;
        gaussian1 = -1.*exp(-1.*((X(:,1)+offset).^2 + (X(:,2)+offset).^2)./100);
        gaussian1 = cat(2,gaussian1,gaussian1);
        
        gaussian2 = exp(-1.*((X(:,1)-offset).^2 + (X(:,2)-offset).^2)./100);
        gaussian2 = cat(2,gaussian2,gaussian2);
        
        grad =  gaussian1.*(offset+X)+gaussian2.*(offset-X);
    case '1gaussian'
        offset = 0;
        gaussian1 = -1.*exp(-1.*((X(:,1)+offset).^2 + (X(:,2)+offset).^2)./100);
        gaussian1 = cat(2,gaussian1,gaussian1);
        grad =  gaussian1.*(offset+X);
    case 'quadraticTrap'
        grad = -2.*X;
    case '2gaussiansOppSign'
        offset = 20;
        gaussian1 = exp(-1.*((X(:,1)+offset).^2 + (X(:,2)+offset).^2)./100);
        gaussian1 = cat(2,gaussian1,gaussian1);
        
        gaussian2 = exp(-1.*((X(:,1)-offset).^2 + (X(:,2)-offset).^2)./100);
        gaussian2 = cat(2,gaussian2,gaussian2);
        grad = gaussian1.*(offset+X)+gaussian2.*(offset-X);
    case 'manyPeaks'
        exponential = exp(-1.*(X(:,2).^2)./4);
        grad = [8.*exponential.*cos(2.*X(:,1)), -2.*X(:,2).*exponential.*sin(2.*X(:,1))];
        grad = grad./10;
    case 'fence'
        a = 0.44;
        inside = (X(:,1)<20).*(X(:,2)<20);
        indices = find(inside);
        grad(indices,:) = (-400/a^2).*[exp((-1/a).*(X(indices,1)-20).^2) exp((-1/a).*(X(indices,2)-20).^2)];
        indices = find(~inside);
        grad(indices,:) = (400/a^2).*[exp((-1/a).*(X(indices,1)-20).^2) exp((-1/a).*(X(indices,2)-20).^2)];
        indicesToZero = (X(:,1)<25).*(X(:,1)>15).*(X(:,2)<25).*(X(:,2)>15);
        grad(find(indicesToZero),:) = 0;
    case 'fenceNoGap'
        a = 0.44;
        inside = (X(:,1)<20).*(X(:,2)<20);
        indices = find(inside);
        grad(indices,:) = (-400/a^2).*[exp((-1/a).*(X(indices,1)-20).^2) exp((-1/a).*(X(indices,2)-20).^2)];
        indices = find(~inside);
        grad(indices,:) = (400/a^2).*[exp((-1/a).*(X(indices,1)-20).^2) exp((-1/a).*(X(indices,2)-20).^2)];
    case 'infiniteFence'
        a = 0.44;
        grad(:,1) = (-400/a^2).*exp((-1/a).*(X(:,1)-20).^2);
    case 'zero'
        return;
end
pert = zeros(numel(X)/2,2);
if ~(strcmp('fence',scalarField))
    for n = 1:numel(c(:,1,1))
        for m = 1:numel(c(1,:,1))
            pert = pert + (1/(n*n*m*m)).*c(n,m,1).*[n.*cos(n.*X(:,1)).*sin(m.*X(:,2)), m.*cos(n.*X(:,1)).*sin(m.*X(:,2))];
            pert = pert + (1/(n*n*m*m)).*c(n,m,2).*[n.*cos(n.*X(:,1)).*cos(m.*X(:,2)), -1.*m.*sin(n.*X(:,1)).*sin(m.*X(:,2))];
            pert = pert + (1/(n*n*m*m)).*c(n,m,3).*[-1.*n.*sin(n.*X(:,1)).*cos(m.*X(:,2)), -1.*m.*cos(n.*X(:,1)).*sin(m.*X(:,2))];
        end
    end
    pert = pert./(numel(c(:,1,1)));
end
%disp(grad)
%quiver(X(:,1),X(:,2),pert(:,1),pert(:,2))
grad = grad + pert;

