function pol = polarization(U,Ndogs)

normedDisplacement = U./vecnorm(U,2,2);
Nsheep = numel(U)/2 - Ndogs;
pol = vecnorm(sum(normedDisplacement(Ndogs+1:end,:)),2,2)/Nsheep;

