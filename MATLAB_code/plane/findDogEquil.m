function V = findDogEquil(CM, L, tar,NumSheepNbhd)

arguments
    CM (1,2) double   % Local Center of Mass
    L double        % attraction radius between two agents
    tar Target      % Target the dogs are going to
    NumSheepNbhd    % # of sheep in the nbhd of the dog
end

% Define constants
R = -1*sqrt(NumSheepNbhd)*L;

% Find equilibrium position (colinear point with CM and dogtarget which
% places CM between V and dog target
V = homeToTarget(tar, CM);
Vnormed = vecnorm(V,2,2);
if ~(Vnormed==0)
    V = V./vecnorm(V,2,2);
end
%scaling taken from observed results
V = V.*R;
V = V + CM;

end