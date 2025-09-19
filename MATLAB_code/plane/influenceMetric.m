function C = influenceMetric(U_T,Tau,t,N,Tmax)    

U_TNorm = vecnorm(U_T,2,2);
U_TNorm(U_TNorm==0) = 1;
U_T = U_T./U_TNorm;
C = zeros(N,N);
for i=1:N
    for j=1:N
        for tPrime = ceil(t-Tau/2):floor(t+Tau/2)
            if (tPrime<1 || tPrime+Tau<1 || tPrime>Tmax || tPrime+Tau>Tmax)
                continue
            end
            vi = U_T(i,:,tPrime);
            vj = U_T(j,:,tPrime+Tau);
            C(i,j) = C(i,j) + dot(vi,vj)/Tau;
        end
    end
end


end