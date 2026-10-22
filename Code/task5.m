clc
clear
close all

%% Task 5 


% Note that while before robustness implementation and verifications worked
% assuming uniform distributions, now we exploit statistical infos about
% our uncertainty to provide insight on the real robustness.

  % task4;
   
  % save('refWorkspace_t4.mat');
   
 clc
 close all

 load('refWorkspace_t4.mat');

%% Complete uncertain model
%  Build a complete model for the quadrotor dynamics using the available 
%  information and the MATLAB Robust Control Toolbox for uncertain system 
%  representation, using:
% • Gaussian densities for the uncertain parameters
% • Uniform density for 𝜓.

%%

% The total system complementary sensitivity transfer function is for
% a MIMO system with 4 inputs and 4 outputs: so it is a 4x4 matrix

% Choose wether to simulate the coupled or uncoupled system:
coupled = 1;

% A. Dyn. performances --------------------------------------------------

% Performance Weights defined in task 1:
OneOverWeightSphi   = 1/W_p_lat; % performance constraints on S
OneOverWeightStheta = 1/W_p_lon;

% Consider the uncertain complete coupled system of task 4.
% The complementary sensitivity functions can be retrieved as:
% I consider the model with the uncertaintes of the plants still lumped
% inside the plant blocks, i.e. without pulling out the deltas, therefore

% we define the imputs:
G_lat_unc.u = 'u_lat';
G_lon_unc.u = 'u_lon';
% tf phi0 -> phi        % uss
F_phi0phi = connect(SumOuterLat,R_phi_tuned, R_p_tuned, G_lat_unc, ... 
            Wc_blk, Delta_c, SumLat, SumLong, ...
            SumOuterLon,R_theta_tuned,R_q_tuned, G_lon_unc, ... 
            'phi_0','phi'); % uss object: uncertain state space (uss) model

% tf theta0 -> theta     % uss
F_theta0theta = connect(SumOuterLat,R_phi_tuned, R_p_tuned, G_lat_unc, ...
            Wc_blk, Delta_c, SumLat, SumLong, ...
            SumOuterLon,R_theta_tuned,R_q_tuned, G_lon_unc, ...
           'theta_0','theta'); 

% That are the diagonal components of the MIMO system's 2x2 F:
F = connect(SumOuterLat,R_phi_tuned, R_p_tuned, G_lat_unc, ... % uss        (G_lat and G_lon are uncertain)
            Wc_blk, Delta_c, SumLat, SumLong, ...
            SumOuterLon,R_theta_tuned,R_q_tuned, G_lon_unc, ...
           {'phi_0','theta_0'},{'phi','theta'});        

% Note that theta_0 is not used in the F_phi0phi and phi_0 is not used in
% F_theta0theta (as shown in the Command Window output).

% Dimensions check: 
% F_phi0phi        1x1 ok
% F_theta0theta    1x1 ok
% F                2x2 ok.

% Check that F_phi0phi is F11 and F_theta0theta is F22:
F_nom = F.NominalValue;  tf(F_nom)            % ss nominale 2x2
F_phi0phi_nom = F_phi0phi.NominalValue;  tf(F_phi0phi_nom) 
F_theta0theta_nom = F_theta0theta.NominalValue; tf(F_theta0theta_nom)
% Ok. and we see that as expected in the nominal case the two dynamics are
% not coupled because the coupling matrix C becomes diagonal, and therefore
% so does F.

% B. Control effort limitation ------------------------------------------

OneOverWeightQphi   = 1/W_q_lat;
OneOverWeightQtheta = 1/W_q_lon;

% tf phi0 -> u_lat
Q_phi0ulat  = connect(SumOuterLat, R_phi_tuned  , R_p_tuned, G_lat_unc, ... % uss (G_lat and G_lon are uncertain)
                        Wc_blk, Delta_c, SumLat, SumLong, ...
                        SumOuterLon, R_theta_tuned, R_q_tuned, G_lon_unc, ...
                       'phi_0','delta_lat'); % uss object

% tf theta0 -> delta_long
Q_theta0ulong = connect(SumOuterLat, R_phi_tuned  , R_p_tuned, G_lat_unc, ... % uss (G_lat and G_lon are uncertain)
                        Wc_blk, Delta_c, SumLat, SumLong, ...
                        SumOuterLon, R_theta_tuned, R_q_tuned, G_lon_unc, ...
                       'theta_0','delta_lon'); % uss object

% Performance verification method

% - We cannot use as performance verification method the evaluation of F peak
% and 0dB crossing frequency because that only happens if F has a
% xi<1/sqrt(2). In all other cases the 2nd order F remains below 0 dB so
% these quantities do not exist. The same stands for an any order F
% function: some F magnitude plots are going to be below 0dB axis, and
% therefore finding the 'peak' and the F-crossing frequency' is not a
% suitable performance verification procedure.

% We use as parameters the same quantities of task 1:
%  - Settling time   
%  - percentage overshoot    
%  - maximum absolute control effort as response to the doublet.



%% Verify via Monte Carlo - estimates the mean values

% Implement Monte Carlo simulation on the complete uncertain feedback 
% system and analyse compliance with the design requirements justifying 
% the choice of the number of MC iterations.


% Choose a k seemingly way smaller than the one we will use
k = 500;


% doublet angle (phi or theta) (variable setpoint of the simulation):
t = linspace(0,6,10^4);
y0 = 0*(t<=1) + 10*(t>1 & t<= 3) - 10*(t>3 & t<= 5) + 0*(t>5); % step [°], doublet setpoint (phi0 or theta0)
y0 = deg2rad(y0); % [rad]

D = Aircraft.Stability;
C = Aircraft.Control;

% sigma = |nominal| * pct/100  (pct = percentual value of sigma wrt nominal value)
sig = @(nom, pct) abs(nom) * pct/100;

% Stability derivatives samples
Yv_samples = normrnd(D.Yv, sig(D.Yv, D.Yv_sigma_pct), [1,k]);
Yp_samples = normrnd(D.Yp, sig(D.Yp, D.Yp_sigma_pct), [1,k]);
Lp_samples = normrnd(D.Lp, sig(D.Lp, D.Lp_sigma_pct), [1,k]);
Xu_samples = normrnd(D.Xu, sig(D.Xu, D.Xu_sigma_pct), [1,k]);
Xq_samples = normrnd(D.Xq, sig(D.Xq, D.Xq_sigma_pct), [1,k]);
Mq_samples = normrnd(D.Mq, sig(D.Mq, D.Mq_sigma_pct), [1,k]);

% Control derivatives samples
Yd_samples = normrnd(C.Yd, sig(C.Yd, C.Yd_sigma_pct), [1,k]);
Ld_samples = normrnd(C.Ld, sig(C.Ld, C.Ld_sigma_pct), [1,k]);
Xd_samples = normrnd(C.Xd, sig(C.Xd, C.Xd_sigma_pct), [1,k]);
Md_samples = normrnd(C.Md, sig(C.Md, C.Md_sigma_pct), [1,k]);

delta_psi_samples = unifrnd(-1, 1, [1,k]);

%  Monte Carlo loop 

% aleatory variables
wc_phi_vec = zeros(1, k);        % crossover frequency
wc_theta_vec = zeros(1, k);      % crossover frequency
xi_phi_vec = zeros(1, k);        % damping
xi_theta_vec = zeros(1, k);      % damping
u_lat_max_vec = zeros(1, k);     % max. control effort in absolute value
u_lon_max_vec = zeros(1, k);     % max. control effort in absolute value
settlTime_phi_vec = zeros(1,k);  % settling time of phi to the step
settlTime_the_vec = zeros(1,k);  % settling time of theta to the step
overshoot_phi_vec = zeros(1, k); % overshoot of phi to the step
overshoot_the_vec = zeros(1,k);  % overshoot of theta to the step


% progressive sample standard deviations:
sd_wc_phi = zeros(1, k);    
sd_wc_theta = zeros(1, k);
sd_xi_phi = zeros(1, k);
sd_xi_theta = zeros(1, k);
sd_u_lat_max = zeros(1, k);
sd_u_lon_max = zeros(1, k);
sd_settlTime_phi = zeros(1,k); 
sd_settlTime_the = zeros(1,k); 
sd_overshoot_phi = zeros(1, k); 
sd_overshoot_the = zeros(1,k); 

% progressive sample means:
mu_wc_phi        = zeros(1, k);
mu_wc_theta      = zeros(1, k);
mu_xi_phi        = zeros(1, k);
mu_xi_theta      = zeros(1, k);
mu_u_lat_max     = zeros(1, k);
mu_u_lon_max     = zeros(1, k);
mu_settlTime_phi = zeros(1,k); 
mu_settlTime_the = zeros(1,k); 
mu_overshoot_phi = zeros(1, k); 
mu_overshoot_the = zeros(1,k);

% magnitude Bode plots of F_phi and F_p
N_om = 500;
omega = logspace(-2, 3, N_om);   % frequency discretization

magFphi = zeros(N_om, k);      % magnitude bode plots of F comp. sensitivities
magFtheta = zeros(N_om, k);
magSphi = zeros(N_om, k);      % magnitude bode plots of S sensitivities
magStheta = zeros(N_om, k);
magQphi = zeros(N_om, k);      % magnitude bode plots of Q control sensitivities
magQtheta = zeros(N_om, k);

for i = 1:k

params_i = struct(...
    'Yv_real',   Yv_samples(i), ...
    'Yp_real',   Yp_samples(i), ...
    'Lp_real',   Lp_samples(i), ...
    'Xu_real',   Xu_samples(i), ...
    'Xq_real',   Xq_samples(i), ...
    'Mq_real',   Mq_samples(i), ...
    'Yd_real',   Yd_samples(i), ...
    'Ld_real',   Ld_samples(i), ...
    'Xd_real',   Xd_samples(i), ...
    'Md_real',   Md_samples(i), ...
    'delta_psi', delta_psi_samples(i)*coupled); 

    % i-th sampled system's tf:
    % usubs: substitutes given values for uncertain elements of uncertain objects
    F_phi0phi_i         = usubs(F_phi0phi        , params_i); % ss (not uss!): not uncertain.
    F_theta0theta_i     = usubs(F_theta0theta    , params_i);
    Q_phi0ulat_i        = usubs(Q_phi0ulat       , params_i);
    Q_theta0ulong_i     = usubs(Q_theta0ulong    , params_i);

    % Compute i-th Lateral and  Longitudinal performances and control effort limitation

    % A. Nominal performance ---------------------------------------------
      % Define a function that, given F, computes wc and xi: see function 
      % compute_wc_xi() at the end of the file. So:
    [wc_phi_vec(i), xi_phi_vec(i), magFphi(:,i), magSphi(:,i)]         = compute_wc_xi(F_phi0phi_i , omega );     % lateral
    [wc_theta_vec(i), xi_theta_vec(i), magFtheta(:,i), magStheta(:,i)] = compute_wc_xi(F_theta0theta_i ,omega );  % longitudinal
    
    stepPhi   = stepinfo(F_phi0phi_i);
    stepTheta = stepinfo(F_theta0theta_i);
    settlTime_phi_vec (i) = stepPhi.SettlingTime;
    settlTime_the_vec (i) = stepTheta.SettlingTime;
    overshoot_phi_vec (i) = stepPhi.Overshoot;
    overshoot_the_vec (i) = stepTheta.Overshoot;

    % B. Control effort limitation ---------------------------------------
      [u_lat_i,tOutLat] = lsim(Q_phi0ulat_i,y0,t);    % 
      [u_lon_i,tOutLon] = lsim(Q_theta0ulong_i,y0,t); % 
       u_lat_max_vec(i) = max(abs(u_lat_i));          %
       u_lon_max_vec(i) = max(abs(u_lon_i));          % 
      [magQphiMAT, ~, ~] = bode(Q_phi0ulat_i, omega); % gives a matrix
       magQphi(:,i) = 20*log10(squeeze(magQphiMAT));  % squeeze, dB
      [magQthetaMAT, ~, ~] = bode(Q_theta0ulong_i, omega);
       magQtheta(:,i) = 20*log10(squeeze(magQthetaMAT));


    % compute sample standard deviations for the aleatory variables after i samples
    sd_wc_phi(i)        = std(wc_phi_vec(1:i));
    sd_wc_theta(i)      = std(wc_theta_vec(1:i));
    sd_xi_phi(i)        = std(xi_phi_vec(1:i));
    sd_xi_theta(i)      = std(xi_theta_vec(1:i));
    sd_u_lat_max(i)     = std(u_lat_max_vec(1:i));
    sd_u_lon_max(i)     = std(u_lon_max_vec(1:i));
    sd_settlTime_phi(i) = std(settlTime_phi_vec(1:i)); 
    sd_settlTime_the(i) = std(settlTime_the_vec(1:i));
    sd_overshoot_phi(i) = std(overshoot_phi_vec(1:i));
    sd_overshoot_the(i) = std(overshoot_the_vec(1:i));

    % compute sample means for the aleatory variables after i samples
    mu_wc_phi(i)        = mean(wc_phi_vec(1:i));
    mu_wc_theta(i)      = mean(wc_theta_vec(1:i));
    mu_xi_phi(i)        = mean(xi_phi_vec(1:i));
    mu_xi_theta(i)      = mean(xi_theta_vec(1:i));
    mu_u_lat_max(i)     = mean(u_lat_max_vec(1:i));
    mu_u_lon_max(i)     = mean(u_lon_max_vec(1:i));
    mu_settlTime_phi(i) = mean(settlTime_phi_vec(1:i));  
    mu_settlTime_the(i) = mean(settlTime_the_vec(1:i));
    mu_overshoot_phi(i) = mean(overshoot_phi_vec(1:i));
    mu_overshoot_the(i) = mean(overshoot_the_vec(1:i));

    fprintf('Sigma it.  %d of %d \n',i,k)

end

sd_data = {sd_wc_phi,        sd_wc_theta, ...
           sd_xi_phi,        sd_xi_theta, ...
           sd_u_lat_max,     sd_u_lon_max,...
           sd_settlTime_phi, sd_settlTime_the, ...
           sd_overshoot_phi, sd_overshoot_the};
%%

MCplot_stdConvergence(sd_data,k)

%%
% !!! Note: to see the gaussianity of the mean of sigma as N->infty you should
% run this MC with the same k for N times and each time evaluate the sample
% sigma with k samples, only at the end, NOT PROGRESSIVE. At the end you
% will have a sigma vector which is N<1 with each sigma computed with k
% samples, and THAT population will have gaussian distribution! But then
% the problem is: how do you choose N? So the same problem you want to
% solve represents itself. If it was solvable by saying ' I just look at
% N s.t. it becomes gaussian enough' then you could just do this directly
% for the mean.
% Here, instead, we run only once and wait for the sigma to converge to a
% value and we use that value as if true sd, even though it is only one of
% the aforementioned samples of a gaussian distribution of sigmas. This
% method that we will follow is the approach also considered by
% Sheldon-Ross. 

%%      

sigma_wc_phi     = sd_wc_phi(end);
sigma_wc_theta   = sd_wc_theta(end);
sigma_xi_phi     = sd_xi_phi(end);
sigma_xi_theta   = sd_xi_theta(end);
sigma_delta_lat  = sd_u_lat_max(end);
sigma_delta_lon  = sd_u_lon_max(end);

% Choice of n for MC simulation -------------------------------------------
% We choose to run a MC for n = k + r iterations (in addition to previous
% ones (note: if you follow Ross's approach the initial k mu values must be
% computed in the previous k iteration, to have less computations, and only
% the remaining r here!)) s.t. "with 90 percent confidence" we assert that the true mean lies
% within 1.96 σ /√n of the observed sample mean. The interval
% is therefore a 95% confidence interval (CI) estimate of μ.
% n is high enough s.t. the sample mean converges.

% Relative error choice ---------------------------------------------------
% Choose a relative error for the performance and control variables:
% we want our final result to be 95% probable of being inside this relative
% error. To respect the orders of magnitude, we consider it as a fraction
% , 1%, of the last sample mean values computed in the previous MC loop (e.g.
% a wc error of 0.1 rad/s). If high level engineering constraint are
% available, then those could be selected instead as errors of wc, xi, etc.
% for which we want 90% confidence. Therefore:

sigmas = [sigma_wc_phi,          sigma_wc_theta, ...
          sigma_xi_phi,           sigma_xi_theta, ...
          sigma_delta_lat,        sigma_delta_lon, ...
          sd_settlTime_phi(end),  sd_settlTime_the(end), ...
          sd_overshoot_phi(end),  sd_overshoot_the(end)];
means  = [mean(wc_phi_vec),        mean(wc_theta_vec), ...
          mean(xi_phi_vec),         mean(xi_theta_vec), ...
          mean(u_lat_max_vec),      mean(u_lon_max_vec), ...
          mean(settlTime_phi_vec),  mean(settlTime_the_vec), ...
          mean(overshoot_phi_vec),  mean(overshoot_the_vec)];
err    = [0.1,0.1,0.01,0.01,0.01,0.01, 0.001,0.001, 0.01,0.01]; 
% choice depends on varible typ. value.

fprintf('\nRequired n for 90%% CI within chosen errors on the mean:\n');
% To satisfy these absolute (not relative) error requirements WITHIN the
% 90% confidence level, each i-th parameter requires a minimum n_i of
% iterations for the convergence, given by the CI semi-width:
% hp: the mean i will compute in the end is a gaussian variable (valid if
% N>> for CLT)
    n = max(max(ceil((1.645 * sigmas ./ err).^2)), k);  % scalar, at least k 

    fprintf(' n >= %d\n', n);
% Therefore we need to simulate at least N simulations, being to satisfy
% the tightiest requirement of confidence between the variables.

% Since we already computed k MC simulation for the mean, then only 
r = n - k;
% are missing. Therefore, once again let's restart the MC simulation loop,
% but without computing the sd (not necessary anymore), only the sample
% means:
%%
% allocate additional r elements:

% Generate ALL additional samples before the loop: append r additional 
% samples to existing vectors:
Yv_samples  = [Yv_samples,  normrnd(D.Yv, sig(D.Yv, D.Yv_sigma_pct), [1,r])];
Yp_samples  = [Yp_samples,  normrnd(D.Yp, sig(D.Yp, D.Yp_sigma_pct), [1,r])];
Lp_samples  = [Lp_samples,  normrnd(D.Lp, sig(D.Lp, D.Lp_sigma_pct), [1,r])];
Xu_samples  = [Xu_samples,  normrnd(D.Xu, sig(D.Xu, D.Xu_sigma_pct), [1,r])];
Xq_samples  = [Xq_samples,  normrnd(D.Xq, sig(D.Xq, D.Xq_sigma_pct), [1,r])];
Mq_samples  = [Mq_samples,  normrnd(D.Mq, sig(D.Mq, D.Mq_sigma_pct), [1,r])];
Yd_samples  = [Yd_samples,  normrnd(C.Yd, sig(C.Yd, C.Yd_sigma_pct), [1,r])];
Ld_samples  = [Ld_samples,  normrnd(C.Ld, sig(C.Ld, C.Ld_sigma_pct), [1,r])];
Xd_samples  = [Xd_samples,  normrnd(C.Xd, sig(C.Xd, C.Xd_sigma_pct), [1,r])];
Md_samples  = [Md_samples,  normrnd(C.Md, sig(C.Md, C.Md_sigma_pct), [1,r])];
delta_psi_samples = [delta_psi_samples, unifrnd(-1, 1,  [1,r])];

% aleatory variables
wc_phi_vec        = [wc_phi_vec,        zeros(1, r)];
wc_theta_vec      = [wc_theta_vec,      zeros(1, r)];
xi_phi_vec        = [xi_phi_vec,        zeros(1, r)];
xi_theta_vec      = [xi_theta_vec,      zeros(1, r)];
u_lat_max_vec     = [u_lat_max_vec, zeros(1, r)];
u_lon_max_vec     = [u_lon_max_vec, zeros(1, r)];
settlTime_phi_vec = [settlTime_phi_vec, zeros(1, r)];
settlTime_the_vec = [settlTime_the_vec, zeros(1, r)];  
overshoot_phi_vec = [overshoot_phi_vec, zeros(1, r)];
overshoot_the_vec = [overshoot_the_vec, zeros(1, r)];  

mu_wc_phi        = [mu_wc_phi,        zeros(1, r)];
mu_wc_theta      = [mu_wc_theta,      zeros(1, r)];
mu_xi_phi        = [mu_xi_phi,        zeros(1, r)];
mu_xi_theta      = [mu_xi_theta,      zeros(1, r)];
mu_u_lat_max     = [mu_u_lat_max, zeros(1, r)];
mu_u_lon_max     = [mu_u_lon_max, zeros(1, r)];
mu_settlTime_phi = [mu_settlTime_phi, zeros(1, r)];
mu_settlTime_the = [mu_settlTime_the, zeros(1, r)];
mu_overshoot_phi = [mu_overshoot_phi, zeros(1, r)];
mu_overshoot_the = [mu_overshoot_the, zeros(1, r)];

magFphi   = [magFphi,   zeros(N_om, r)];
magFtheta = [magFtheta, zeros(N_om, r)];
magSphi   = [magSphi,   zeros(N_om, r)];
magStheta = [magStheta, zeros(N_om, r)];
magQphi   = [magQphi,   zeros(N_om, r)];
magQtheta = [magQtheta, zeros(N_om, r)];


for i = 1:r

 
params_i = struct(...
    'Yv_real',   Yv_samples(k + i), ...
    'Yp_real',   Yp_samples(k + i), ...
    'Lp_real',   Lp_samples(k + i), ...
    'Xu_real',   Xu_samples(k + i), ...
    'Xq_real',   Xq_samples(k + i), ...
    'Mq_real',   Mq_samples(k + i), ...
    'Yd_real',   Yd_samples(k + i), ...
    'Ld_real',   Ld_samples(k + i), ...
    'Xd_real',   Xd_samples(k + i), ...
    'Md_real',   Md_samples(k + i), ...
    'delta_psi', delta_psi_samples(k + i)*coupled);

    % i-th sampled system's tf:
    % usubs: substitutes given values for uncertain elements of uncertain objects
    F_phi0phi_i         = usubs(F_phi0phi        , params_i); % ss (not uss!): not uncertain.
    F_theta0theta_i     = usubs(F_theta0theta    , params_i);
    Q_phi0ulat_i        = usubs(Q_phi0ulat       , params_i);
    Q_theta0ulong_i     = usubs(Q_theta0ulong    , params_i);

    % Compute i-th Lateral and  Longitudinal performances and control effort limitation

    % A. Nominal performance ---------------------------------------------
      % Define a function that given F, computes wc and xi: see function 
      % compute_wc_xi() at the end of the file. So:
    [wc_phi_vec(k + i), xi_phi_vec(k + i), magFphi(:,k+i), magSphi(:,k+i)]     = compute_wc_xi(F_phi0phi_i, omega );         % lateral
    [wc_theta_vec(k + i), xi_theta_vec(k + i), magFtheta(:,k+i), magStheta(:,k+i)] = compute_wc_xi(F_theta0theta_i, omega ); % longitudinal

    stepPhi   = stepinfo(F_phi0phi_i);
    stepTheta = stepinfo(F_theta0theta_i);
    settlTime_phi_vec (k + i) = stepPhi.SettlingTime;
    settlTime_the_vec (k + i) = stepTheta.SettlingTime;
    overshoot_phi_vec (k + i) = stepPhi.Overshoot;
    overshoot_the_vec (k + i) = stepTheta.Overshoot;

    % B. Control effort limitation ---------------------------------------
      [u_lat_i,tOutLat] = lsim(Q_phi0ulat_i,y0,t);    % simulate roll (phi) doublet response
      [u_lon_i,tOutLon] = lsim(Q_theta0ulong_i,y0,t); % simulate pitch (theta) doublet response
       u_lat_max_vec(k + i) = max(abs(u_lat_i));      % max.(abs.) control action to respond to the roll  doublet
       u_lon_max_vec(k + i) = max(abs(u_lon_i));      % max.(abs.) control action to respond to the pitch doublet
      [magQphiMAT, ~, ~]  = bode(Q_phi0ulat_i, omega);
       magQphi(:,k+i) = 20*log10(squeeze(magQphiMAT));
      [magQthetaMAT, ~, ~] = bode(Q_theta0ulong_i, omega);
       magQtheta(:,k+i) = 20*log10(squeeze(magQthetaMAT));

    % compute sample means for the aleatory variables after i samples
    mu_wc_phi(k + i)        = mean(wc_phi_vec(1:k + i));
    mu_wc_theta(k + i)      = mean(wc_theta_vec(1:k + i));
    mu_xi_phi(k + i)        = mean(xi_phi_vec(1:k + i));
    mu_xi_theta(k + i)      = mean(xi_theta_vec(1:k + i));
    mu_u_lat_max(k + i)     = mean(u_lat_max_vec(1:k + i));
    mu_u_lon_max(k + i)     = mean(u_lon_max_vec(1:k + i));
    mu_settlTime_phi(k + i) = mean(settlTime_phi_vec(1:k + i));
    mu_settlTime_the(k + i) = mean(settlTime_the_vec(1:k + i));
    mu_overshoot_phi(k + i) = mean(overshoot_phi_vec(1:k + i));
    mu_overshoot_the(k + i) = mean(overshoot_the_vec(1:k + i));

    fprintf('Mu it.  %d of %d \n',i,r)

end

mu_data = {mu_wc_phi,        mu_wc_theta, ...
           mu_xi_phi,        mu_xi_theta, ...
           mu_u_lat_max,     mu_u_lon_max, ...
           mu_settlTime_phi, mu_settlTime_the, ...
           mu_overshoot_phi, mu_overshoot_the};

%%
% !!!: Also store all the step responses values and print the plot of all
% step responses like they do with the Robust Control Toolbox!

% !!! : Mostra anche che le Q e la F (la S actually, = 1-F easily computable)
%  stanno sotto i pesi scegli nel nominal design.

% N.B. you cannot compute with a for loop like this one the uncertainty 
% regions in the polar plot because due to the gaussianness of the variable 
% technically they would be infinite, but 'improbable' exponentially from 
% the nominal polar plot points. You have to consider e.g. constant
% distribution within 3 sigma (as you also defined the ureal variables that
% should actualy be gaussian). To do it you should generate them as
% uniformly distributed in the MC and then you select some omegas at which
% you want to build the uncertainty regions on the nominal polar plot
% curve.

MCplot_meanConvergence(mu_data,n)

%% Bode plots MC envelopes VS Nominal

% Performances ------------------------------------------------------------
% Pre-compute weight magnitudes
[magWS_phi,   ~, ~] = bode(OneOverWeightSphi,   omega);
[magWS_theta, ~, ~] = bode(OneOverWeightStheta, omega);
[magWQ_phi,   ~, ~] = bode(OneOverWeightQphi,   omega);
[magWQ_theta, ~, ~] = bode(OneOverWeightQtheta, omega);
magWS_phi   = 20*log10(squeeze(magWS_phi));
magWS_theta = 20*log10(squeeze(magWS_theta));
magWQ_phi   = 20*log10(squeeze(magWQ_phi));
magWQ_theta = 20*log10(squeeze(magWQ_theta));

%  - - - - - F_phi - - - - -
S_phi0phi_nom = 1 - F_phi0phi_nom;
[magF_phi_nom, ~, ~] = bode(F_phi0phi_nom, omega);
[magS_phi_nom, ~, ~] = bode(S_phi0phi_nom, omega);
magF_phi_nom = 20*log10(squeeze(magF_phi_nom));
magS_phi_nom = 20*log10(squeeze(magS_phi_nom));

figure;
for j = 1:n
    semilogx(omega, magFphi(:,j), 'Color', [0.27 0.48 0.78 0.15], 'LineWidth', 0.1); hold on
    semilogx(omega, magSphi(:,j), 'Color', [0.90 0.22 0.27 0.15], 'LineWidth', 0.1);
end
hFphi  = semilogx(omega, magF_phi_nom, 'Color', [0.27 0.48 0.78], 'LineWidth', 0.1);
hSphi  = semilogx(omega, magS_phi_nom, 'Color', [0.90 0.22 0.27], 'LineWidth', 0.1);
hWSphi = semilogx(omega, magWS_phi,    'r--', 'LineWidth', 0.1);   % weight bound
yline(0, 'k--', 'LineWidth', 1);
grid on
legend([hFphi hSphi hWSphi], {'F_\phi (MC + nominal)', 'S_\phi (MC + nominal)', '1/W_{P,\phi} bound'}, 'Location', 'southwest')
xlabel('\omega [rad/s]'); ylabel('Magnitude [dB]')
title('Bode Magnitude: F_\phi and S_\phi — Monte Carlo envelope')

%  - - - - - F_theta - - - - -
F_theta_nom = F_theta0theta.NominalValue;
S_theta_nom = 1 - F_theta_nom;
[magF_theta_nom, ~, ~] = bode(F_theta_nom, omega);
[magS_theta_nom, ~, ~] = bode(S_theta_nom, omega);
magF_theta_nom = 20*log10(squeeze(magF_theta_nom));
magS_theta_nom = 20*log10(squeeze(magS_theta_nom));

figure;
for j = 1:n
    semilogx(omega, magFtheta(:,j), 'Color', [0.27 0.48 0.78 0.15], 'LineWidth', 0.1); hold on
    semilogx(omega, magStheta(:,j), 'Color', [0.90 0.22 0.27 0.15], 'LineWidth', 0.1);
end
hFtheta  = semilogx(omega, magF_theta_nom,  'Color', [0.27 0.48 0.78], 'LineWidth', 0.1);
hStheta  = semilogx(omega, magS_theta_nom,  'Color', [0.90 0.22 0.27], 'LineWidth', 0.1);
hWStheta = semilogx(omega, magWS_theta,     'r--', 'LineWidth', 0.1);
yline(0, 'k--', 'LineWidth', 1);
grid on
legend([hFtheta hStheta hWStheta], {'F_\theta (MC + nominal)', 'S_\theta (MC + nominal)', '1/W_{P,\theta} bound'}, 'Location', 'southwest')
xlabel('\omega [rad/s]'); ylabel('Magnitude [dB]')
title('Bode Magnitude: F_\theta and S_\theta — Monte Carlo envelope')

% Control effort limitation -----------------------------------------------
%  - - - - - Q_phi - - - - -
Qphi_nom = Q_phi0ulat.NominalValue;
[magQphi_nom, ~, ~] = bode(Qphi_nom, omega);
magQphi_nom = 20*log10(squeeze(magQphi_nom));

figure;
for j = 1:n
    semilogx(omega, magQphi(:,j), 'Color', [0.27 0.48 0.78 0.15], 'LineWidth', 0.1); hold on
end
hQphi  = semilogx(omega, magQphi_nom, 'Color', [0.27 0.48 0.78], 'LineWidth', 0.1);
hWQphi = semilogx(omega, magWQ_phi,   'r--', 'LineWidth', 0.1);
yline(0, 'k--', 'LineWidth', 1);
grid on
legend([hQphi hWQphi], {'Q_\phi (MC + nominal)', '1/W_{Q,\phi} bound'}, 'Location', 'southwest')
xlabel('\omega [rad/s]'); ylabel('Magnitude [dB]')
title('Bode Magnitude: Q_\phi — Monte Carlo envelope')

%  - - - - - Q_theta - - - - -
Qtheta_nom = Q_theta0ulong.NominalValue;
[magQtheta_nom, ~, ~] = bode(Qtheta_nom, omega);
magQtheta_nom = 20*log10(squeeze(magQtheta_nom));

figure;
for j = 1:n
    semilogx(omega, magQtheta(:,j), 'Color', [0.90 0.22 0.27 0.15], 'LineWidth', 0.1); hold on
end
hQtheta  = semilogx(omega, magQtheta_nom, 'Color', [0.90 0.22 0.27], 'LineWidth', 0.1);
hWQtheta = semilogx(omega, magWQ_theta,   'r--', 'LineWidth', 0.1);
yline(0, 'k--', 'LineWidth', 1);
grid on
legend([hQtheta hWQtheta], {'Q_\theta (MC + nominal)', '1/W_{Q,\theta} bound'}, 'Location', 'southwest')
xlabel('\omega [rad/s]'); ylabel('Magnitude [dB]')
title('Bode Magnitude: Q_\theta — Monte Carlo envelope')

% It is possible that some realizations do not respect the requirements because 
% for the gaussian variables the robustness (RP, RS) was implemented s.t.
% the performances are verified for all combinations of variables each
% being able to vary within their 3sigma interval (ureal models as if
% uniformly distributed), so if here we generate the variables within a
% 3sigma variation only (as usample() would), we would get all realizations
% respecting the requirements; but since we here generate the gaussian
% variables as actually gaussian, some of them are generated outside the
% 3sigma threshold (the probability of this happening is < (1-99.73)% =
% 0.27%) FOR each variable. Although note that the probability of not
% verifiying the requirements is not of 0.27%, because the non-fulfillment
% of the requirement could happen also for a case in which all the
% slightly-off realizations of the 11 variables, even if within their own
% 3sigma threshold, combine together to make the system change quite a lot.
% So it is possible (and with N = infty iterations it is sure), and expected,
% that some realizations will not satisfy the requirements, due to the
% fact that the modelling considered only uniform distributions (3sigma
% is considered the engineering 0 for probability distribution functions),
% while the real system has gaussian uncertaintes, that possibly can have
% the uncertain values be any real value in fact (gauss bell doesn't ever
% reach 0).
%  
 

%% Functions

% Define a function that given F, computes wc and xi:

function [wc, xi, magF, magS] = compute_wc_xi(F, omega)

    % --- Recover loop gain from F ---
    L = F / (1 - F); 

    % --- Bode of L ---
    [mag, phase, wout] = bode(L, omega);
    mag   = 20*log10(squeeze(mag));   % [dB]
    phase = squeeze(phase);           % [deg] 
    wout  = squeeze(wout);                  

    % --- A.1 Crossover frequency: |L|_dB = 0 dB ---
    wc = interp1(mag, wout, 0);

    % --- A.2 Crossover phase and phase margin ---
    % Since F is not necessary of order 2, the definition of xi is
    % through its 2nd order representation:
    phi_c = interp1(wout, phase, wc); 
    phi_m = 180 - abs(phi_c);         % phase margin [deg]
    xi = sin(deg2rad(phi_m) / 2);     
    % alternatively: xi ≈ phi_m / 100 
     
    % F bode mag plot
    [magF, ~, ~] = bode(F, omega);
    magF = 20*log10(squeeze(magF));

    % S bode mag plot
    S = 1 - F;
    [magS, ~, ~] = bode(S, omega);
    magS = 20*log10(squeeze(magS));

end

%% Plots

function MCplot_stdConvergence(sd_data, k)

figure('Position', [100, 100, 1400, 700]);

titles = {'\omega_c^\phi [rad/s]',   '\omega_c^\theta [rad/s]', ...
          '\xi_\phi [-]',            '\xi_\theta [-]', ...
          '\delta_{lat}^{max} [rad]', '\delta_{lon}^{max} [rad]', ...
          'T_{set}^\phi [s]',        'T_{set}^\theta [s]', ...
          'OS_\phi [%]',             'OS_\theta [%]'};

colors = {'#E63946','#457B9D','#2A9D8F','#E9C46A','#F4A261','#A8DADC', ...
          '#6A0572','#1D7A8A','#E07B39','#3D9970'};

for j = 1:10
    subplot(2, 5, j);
    plot(1:k, sd_data{j}, 'Color', colors{j}, 'LineWidth', 1.5);
    hold on
    yline(sd_data{j}(end), 'k--', 'LineWidth', 1.2, ...
          'Label', sprintf('\\sigma_{\\infty}=%.4f', sd_data{j}(end)), ...
          'LabelOrientation', 'horizontal');
    title(titles{j}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('N iterations');
    ylabel('Progressive \sigma');
    grid on; box on;
end

sgtitle(sprintf('Progressive Std Dev — CLT Convergence  (N = %d)', k), ...
        'FontSize', 13, 'FontWeight', 'bold');
end


function MCplot_meanConvergence(mu_data, n)

figure('Position', [100, 100, 1400, 700]);

titles = {'\omega_c^\phi [rad/s]',   '\omega_c^\theta [rad/s]', ...
          '\xi_\phi [-]',            '\xi_\theta [-]', ...
          '\delta_{lat}^{max} [rad]', '\delta_{lon}^{max} [rad]', ...
          'T_{set}^\phi [s]',        'T_{set}^\theta [s]', ...
          'OS_\phi [%]',             'OS_\theta [%]'};

colors = {'#E63946','#457B9D','#2A9D8F','#E9C46A','#F4A261','#A8DADC', ...
          '#6A0572','#1D7A8A','#E07B39','#3D9970'};

for j = 1:10
    subplot(2, 5, j);
    plot(1:n, mu_data{j}, 'Color', colors{j}, 'LineWidth', 1.5);
    hold on
    yline(mu_data{j}(end), 'k--', 'LineWidth', 1.2, ...
          'Label', sprintf('\\mu_{\\infty}=%.4f', mu_data{j}(end)), ...
          'LabelOrientation', 'horizontal');
    title(titles{j}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('N iterations');
    ylabel('Progressive \mu');
    grid on; box on;
end

sgtitle(sprintf('Progressive Means — CLT Convergence  (N = %d)', n), ...
        'FontSize', 13, 'FontWeight', 'bold');
end


%% Plot distributions and CDF (Cumulative Distribution Function)

vars = {wc_phi_vec,       wc_theta_vec, ...
        xi_phi_vec,        xi_theta_vec, ...
        u_lat_max_vec,     u_lon_max_vec, ...
        settlTime_phi_vec, settlTime_the_vec, ...
        overshoot_phi_vec, overshoot_the_vec};

labels = {'\omega_c^\phi  [rad/s]', '\omega_c^\theta  [rad/s]', ...
          '\xi_\phi  [-]',          '\xi_\theta  [-]', ...
          '\delta_{lat}^{max}  [rad]', '\delta_{lon}^{max}  [rad]', ...
          'T_{set}^\phi  [s]',      'T_{set}^\theta  [s]', ...
          'OS_\phi  [%]',           'OS_\theta  [%]'};

colors = {'#E63946','#457B9D','#2A9D8F','#E9C46A','#F4A261', ...
          '#A8DADC','#6A0572','#1D7A8A','#E07B39','#3D9970'};

% ── Figure 1: Histograms ─────────────────────────────────────────────────
figure('Position',[50 50 1400 600]);
sgtitle(sprintf('Monte Carlo (n = %d) — Relative Frequency Histograms', n), ...
        'FontSize',13,'FontWeight','bold')

for j = 1:10
    subplot(2,5,j)
    x = vars{j};
    histogram(x, 35, ...
              'Normalization', 'probability', ...
              'FaceColor',     colors{j}, ...
              'EdgeColor',     'white', ...
              'FaceAlpha',     0.85)
    hold on
    mu_j = mean(x);  s_j = std(x);
    xline(mu_j,       'k-',  'LineWidth',1.8, ...
          'Label', sprintf('\\mu=%.3g', mu_j), ...
          'LabelOrientation','horizontal','LabelVerticalAlignment','bottom')
    xline(mu_j + s_j, 'k--', 'LineWidth',1.0)
    xline(mu_j - s_j, 'k--', 'LineWidth',1.0)
    xlabel(labels{j},'FontSize',9)
    ylabel('Relative freq.','FontSize',9)
    title(labels{j},'FontSize',10,'FontWeight','bold')
    grid on; box on
end

% ── Figure 2: Empirical CDFs ─────────────────────────────────────────────
figure('Position',[50 50 1400 600]);
sgtitle(sprintf('Monte Carlo (n = %d) — Empirical CDF', n), ...
        'FontSize',13,'FontWeight','bold')

for j = 1:10
    subplot(2,5,j)
    x     = sort(vars{j});
    p_cdf = (1:n) / n;
    stairs(x, p_cdf, 'Color', colors{j}, 'LineWidth', 2)
    hold on
    mu_j = mean(vars{j});
    xline(mu_j, 'k-', 'LineWidth',1.5, ...
          'Label', sprintf('\\mu=%.3g', mu_j), ...
          'LabelOrientation','horizontal')
    xlabel(labels{j},'FontSize',9)
    ylabel('CDF','FontSize',9)
    ylim([0 1]); grid on; box on
    title(labels{j},'FontSize',10,'FontWeight','bold')
end

% ── Summary table ────────────────────────────────────────────────────────
fprintf('\n%-22s  %10s  %10s\n','Variable','Mean','Std')
fprintf('%s\n', repmat('-',1,46))
var_names = {'wc_phi','wc_theta','xi_phi','xi_theta', ...
             'u_lat_max','u_lon_max','T_set_phi','T_set_theta', ...
             'OS_phi','OS_theta'};
for j = 1:10
    fprintf('%-22s  %10.4g  %10.4g\n', var_names{j}, mean(vars{j}), std(vars{j}))
end

% Comments: non gaussianity is correct, indeed CLT only constrains the
% samples means to be gaussian. The peak in 0 for the percentage overshoot 
% (OS) tells us the many realizations do not have overshoot at all (OS=0)

% Comparison Coupled vs uncoupled:
% Without coupling, we never have overshoot of phi in any realization.
% Adding coupling, the OS increase, due to coupling evidently. Indeed,
% also the peaks of the bode plots are affected. Also, the
% distributions of some variables become less 'gaussian'.

%% Print CI (Confidence intervals)

%% Print CI (Confidence intervals)
% Compute CI using new sigma, more precise: confidence won't be 90%
% exactly; for some variable it will be more because we used a N more than
% the necessary one!

fprintf('\n%-22s  %10s  %10s  %10s  %10s  %10s\n', ...
        'Variable','Mean','Sigma_N','CI_low','CI_high','Confidence [%]')
fprintf('%s\n', repmat('-',1,76))

sigmas_final = [std(wc_phi_vec),       std(wc_theta_vec), ...
                std(xi_phi_vec),        std(xi_theta_vec), ...
                std(u_lat_max_vec),     std(u_lon_max_vec), ...
                std(settlTime_phi_vec), std(settlTime_the_vec), ...
                std(overshoot_phi_vec), std(overshoot_the_vec)];

for j = 1:10
    mu_j    = mean(vars{j});
    sigma_j = sigmas_final(j);
    e_j     = err(j);

    % quantile z_{1-alpha/2} = e * sqrt(n)/sigma_N
    quantile_j = e_j * sqrt(n) / sigma_j;

    % alpha = 1 - Phi(quantile), two-sided: confidence = 1 - 2*alpha
    alpha_j      = 2*(1 - normcdf(quantile_j));
    confidence_j = (1 - alpha_j) * 100; % defined as 1 - alpha

    CI_low  = mu_j - e_j;
    CI_high = mu_j + e_j;

    fprintf('%-22s  %10.4g  %10.4g  %10.4g  %10.4g  %10.2f%%\n', ...
            var_names{j}, mu_j, sigma_j, CI_low, CI_high, confidence_j)
end