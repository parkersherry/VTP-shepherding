function y = transition(x,spec)

arguments
    x (1,1) double {mustBeNonnegative}
    spec (1,:) char
end

if (x >= 1 && ~strcmp(spec,'linear2'))
    y = 0;
else
    switch spec
        case 'indicator'
            y = 1;
        case 'sin'
            y = sin(pi*x);
        case 'mollifier'
            y = exp(1-1/(1-x^2));
        case 'cubic'
            y = 2*x^3 - 3*x^2 + 1;
        case 'linear'
            y = 1-x;
        case 'linear2'
            y = x^2;
        case 'expReciprocal'
            y = exp(-1/(1-x)) / (exp(-1/x)+exp(-1/(1-x)));
        case 'attractionExp'
            syms z;
            assume(z==x);
            y = piecewise((0 < z) && (z < 0.5), exp(-1/(1-2*z)) / (exp(-1/(2*z))+exp(-1/(1-2*z))), (0.5 < z) && (z < 1), exp(-1/(1-(2*z-1))) / (exp(-1/(2*z-1))+exp(-1/(1-(2*z-1))))-1);
            subs(y,z,x);
        case 'attractionSub'
            y = (1+1)*(exp(-1/(1-x)) / (exp(-1/x)+exp(-1/(1-x))))-1;
        case 'attractionSubAlign'
             y = 1- abs(2*(exp(-1/(1-x)) / (exp(-1/x)+exp(-1/(1-x))))-1);
        case 'dogExpReciprocal'
            y = 3*exp(-1/(1-x)) / (exp(-1/x)+exp(-1/(1-x)));
    end
end