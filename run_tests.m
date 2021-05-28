% This scripts prepares the experimental setup and launches the
% experiments. It is called by the script main.m
% The results are saved in save_name.mat


% create m and N vectors
mm = round(logspace(log10(m_min),log10(m_max),Nmm));  % equally spaced on log scale
NN = round(linspace(Nmin,Nmax,NumbN)); 

% create the statistical model
stat_mod.mu = mu;
stat_mod.k_Sx = k_Sx;
stat_mod.c_k = c_k;
stat_mod.sigma = sigma;
stat_mod.error_model = error_model;

% Store the excess risk for each resolution, sample size and repetition
Gen_noT_par = zeros(length(NN),Nmm,K);
Gen_T_par = zeros(length(NN),Nmm,K);

% Run main code
for iN = 1:length(NN)
    parfor kk = 1:K
        fprintf('N = %d; k = %d/%d \n',NN(iN),kk,K)
        [GU,GS] = generalization_compare(NN(iN),stat_mod,mm);
        Gen_U_par(iN,:,kk) = GU;
        Gen_S_par(iN,:,kk) = GS;
    end
end

% Compute averages over the realization
Gen_U = mean(Gen_U_par,3);
Gen_S = mean(Gen_S_par,3);

% Save
save(save_name,'Gen_U','Gen_S','mm','NN','Gen_U_par','Gen_S_par','stat_mod')


