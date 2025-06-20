function [X,Y,Z] = getMeshScalarField(axes,scalarField)
[X,Y] = meshgrid(axes,axes);
switch scalarField
    case '2gaussiansPositive'
        offset = 20;
        gaussian1 = exp(-1.*((X+offset).^2 + (Y+offset).^2)./100);        
        gaussian2 = exp(-1.*((X-offset).^2 + (Y-offset).^2)./100);
        
        Z =  gaussian1+gaussian2;
    case '1gaussian'
        offset = 0;
        gaussian1 = exp(-1.*((X+offset).^2 + (Y+offset).^2)./100);
        Z =  gaussian1;
    case 'quadraticTrap'
        Z = -1.*(X.^2+Y.^2);
    case '2gaussiansOppSign'
        offset = 20;
        gaussian1 = exp(-1.*((X+offset).^2 + (Y+offset).^2)./100);
        
        gaussian2 = exp(-1.*((X-offset).^2 + (Y-offset).^2)./100);
        Z = gaussian1-gaussian2;
    case 'manyPeaks'
        exponential = exp(-1.*(Y.^2)./4);
        Z = exponential.*sin(2.*X)./10;
    case 'fence'
        Z = zeros(size(X));
        Z(X>20) = 6;
        Z(Y>20) = 6;

        Z(X>21) = 0;
        Z(Y>21) = 0;
        inside = inpolygon(X(:),Y(:),[15;25;25;15],[15;15;20;20]);
        Z(find(inside==1)) = 0;
end