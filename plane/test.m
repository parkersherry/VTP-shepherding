
N = [75 105 135 165 195 225 255 285 315 345];
nu = [0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25 8.75 9.25 9.75];


for nuIndex = 1:20
    posSeed = randi(100000);
    angSeed = randi(100000);
    fnameDogless = '/Users/mikey/Summer 2025/plane/OvernightM26/NoDog/ND_N_250_nu_'+string(nu(nuIndex)) + '_run_' + string(nuIndex) + '.mat';

    serialSandbox(fnameDogless,250,3.3,nu(nuIndex),1000,posSeed,angSeed,0)

    for NIndex = 1:10

        for seedIndex = 1:5
            posSeed = randi(100000);
            angSeed = randi(100000);
            fnameDog = '/Users/mikey/Summer 2025/plane/OvernightM26/Dog/D_N_'+string(N(NIndex)) + '_nu_'+string(nu(nuIndex)) + '_run_' + string(seedIndex) + '.mat';
            serialSandbox(fnameDog,N(NIndex),3.3,nu(nuIndex),1000,posSeed,angSeed,1)

        end


    end
    disp(nuIndex);


end

