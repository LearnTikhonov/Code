function [GeneralizationU,GeneralizationS] = generalization_compare(Np,stat_mod,mm)


% Import the statistical model
mu_fun           = stat_mod.mu;
L_filt_phys_size = stat_mod.c_k;
L_filt_fun       = stat_mod.k_Sx;
sigma            = stat_mod.sigma;
error_model      = stat_mod.error_model;

% Setup
phys_L = 1; % size of the physical domain
phys_x = linspace(0,phys_L,Np); % pixel subdivision
Nmm = length(mm); 

% Random variable definitions

% Definition of x: Gaussian variable with mean and whose covariance is a
% convolution operator L
mu_vec = mu_fun(phys_x)';                 % discretize the mean
L_filt_size = Np*L_filt_phys_size/phys_L; % discretize the kernel support
nu = floor(L_filt_size/2); L_filt_size = 2*nu+1; % use odd-sized filter
L_filt_x = linspace(-L_filt_phys_size,L_filt_phys_size,L_filt_size);
L_filt = L_filt_fun(L_filt_x)*sqrt(phys_L/Np); % discretize the kernel
Sx_half = @(x) circshift(cconv(x,L_filt,Np),-nu); % implements convolution on the torus
Sx_fun = @(x) Sx_half(Sx_half(x));

% Definition of epsilon: white noise process
switch error_model
    case 3 % white noise with respect to the wavelet components
        % Wav  = @(x) haart_vec(x);
        WavT = @(w) haart_vecT(w);
        Se_half = @(x) sigma*sqrt(Np)*WavT(x);
        Se_fun = @(x) sigma^2*Np*x;  % we exploit that Wav(WavT(x))=x
        Se_mat = sigma^2*Np*eye(Np);
        Se_inv_mat = 1/(sigma^2*Np)*eye(Np); 
    otherwise % white noise with respect to the pixel basis
        Se_half = @(x) sigma*sqrt(Np)*x;
        Se_fun = @(x) sigma^2*Np*x;
        Se_mat = sigma^2*Np*eye(Np);
        Se_inv_mat = 1/(sigma^2*Np)*eye(Np);
end

% Forward operator (suited for any diagonal operator)
diagonal_A = ones(Np,1);   % denoising operator
A = @(x) diagonal_A.*x;
At = @(x) diagonal_A.*x; % transpose
A_mat = diag(diagonal_A);

%%% BEST PARAMETER
% MSE computation: tr( W(A Sx A* + Se)W* - 2 Sx A*W* + Sx ) + ||(WA-I)mu + b||^2
%   best parameter: B^2 = Sx, h = mu
C_best = @(y) A(Sx_fun(At(y))) + Se_fun(y);         % auxiliary operator
W_best = @(y) Sx_fun(At(pcg_noWrite(C_best,y)));    % W = Sx A*(A Sx A* + Se)^{-1}
Wt_best = @(y) pcg_noWrite(C_best,A(Sx_fun(y)));    % adjoint of W
% b_best = mu_vec - W_best(A(mu_vec));      % not needed
OP_best = @(y) W_best(C_best(Wt_best(y))) - Sx_fun(At(Wt_best(y))) - W_best(A(Sx_fun(y))) + Sx_fun(y);
MSE_best = trace_op(OP_best,Np)/Np; %+ norm(W_best(A(mu_vec))-mu_vec+b_best)/sqrt(Np);
fprintf('MSE_best = %.6f \n',MSE_best);


%%% TRAINING

GeneralizationU = zeros(1,Nmm);
GeneralizationS = zeros(1,Nmm);

for i = 1:Nmm
    
    % Training sample
    m = mm(i); % sample size
    samp_e = zeros(Np,m);
    samp_x = zeros(Np,m);
    for im = 1:m
        switch error_model
            case 1      % Gaussian noise
                nue = randn(Np,1);
            otherwise   % uniform noise
                nue = (rand(Np,1)-0.5)*sqrt(12); 
        end
        nux = randn(Np,1);
        samp_e(:,im) = Se_half(nue);
        samp_x(:,im) = mu_vec + Sx_half(nux);
    end
    
    %%%%%%%%%%%%%%%%
    % Unsupervised %
    %%%%%%%%%%%%%%%%
    
    % Sample approximation of mu and Sx
    mu_emp = mean(samp_x,2);
    Sx_emp = 1/m*(samp_x-mu_emp*ones(1,m))*(samp_x-mu_emp*ones(1,m))';

    % MSE computation: tr( W(A Sx A* + Se)W* - 2 Sx A*W* + Sx ) + ||(WA-I)mu + b||^2
    %   best parameter: B^2 = Sx_emp, h = mu_emp
    C_u = @(y) A(Sx_emp*(At(y))) + Se_fun(y);
    W_u = @(y) Sx_emp*(At(pcg_noWrite(C_u,y)));
    Wt_u = @(y) pcg_noWrite(C_u,A(Sx_emp*y));
    b_u = mu_emp - W_u(A(mu_emp));
    OP_u = @(y) W_u(C_best(Wt_u(y))) - Sx_fun(At(Wt_u(y))) - W_u(A(Sx_fun(y))) + Sx_fun(y);
    MSE = trace_op(OP_u,Np)/Np + norm(W_u(A(mu_vec))-mu_vec + b_u)/sqrt(Np);

    fprintf('Sample size: m = %d; MSE = %.6f \n',mm(i),MSE);
    GeneralizationU(i) = abs(MSE-MSE_best);    

    
    %%%%%%%%%%%%%%
    % Supervised %
    %%%%%%%%%%%%%%
    
    % Approximate relevant statistical quantities
    eps_emp = mean(samp_e,2);
    Se_emp = 1/m*(samp_e-eps_emp*ones(1,m))*(samp_e-eps_emp*ones(1,m))';
    Sxe_emp = 1/m*(samp_x-mu_emp*ones(1,m))*(samp_e-eps_emp*ones(1,m))'; 
    % Get a symmetric M
    M_star = (Sx_emp*A_mat' + Sxe_emp)*(inv(A_mat'*Se_inv_mat*Se_emp + A_mat'*Se_inv_mat*Sxe_emp'*A_mat'));
    M_star = (M_star+M_star')/2; % take only the symmetric part
    W_star = M_star*A_mat'*(inv(A_mat'*M_star*A_mat + Se_mat));
    b_star = (eye(Np) - W_star*A_mat)*mu_emp - W_star*eps_emp;
   
    % MSE computation: tr( W(A Sx A* + Se)W* - 2 Sx A*W* + Sx ) + ||(WA-I)mu + b||^2
    W_s = @(y) W_star*y;
    Wt_s = @(y) W_star'*y;
    b_s = b_star;
    OP_s = @(y) W_s(C_best(Wt_s(y))) - Sx_fun(At(Wt_s(y))) - W_s(A(Sx_fun(y))) + Sx_fun(y);
    MSE = trace_op(OP_s,Np)/Np + norm(W_s(A(mu_vec))-mu_vec + b_s)/sqrt(Np);
    
    fprintf('  Sample size: m = %d; MSE = %.6f \n',mm(i),MSE);
    GeneralizationS(i) = abs(MSE-MSE_best);
      
end



end
