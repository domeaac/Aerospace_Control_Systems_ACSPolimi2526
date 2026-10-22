clc
clear
close all


%% READ ME FIRST

% In all 6 tasks before this one (1 2 3 lat and long) put clc clear close
% all as commented out.

%% 

    task3_lonGIO; % runs also task1_lon and task2_lon
    
    task3_latGIO; % runs also task1_lat and task2_lat
    
    save('refWorkspace_t1t2t3_latlon.mat');

 clc
 close all

 % load('refWorkspace_t1t2t3_latlon.mat');

%%

s = tf('s');

%% Task 4: Robustness to coupling

% Verify if the obtained control system satisfies the RS condition when 
% taking into account the longitudinal/lateral coupling given by

phi_min = -deg2rad(15);
phi_max =  deg2rad(15);
% OBS.: the value has also a sign that is uncertain. This can cause bvery
% different behaviour in the different realization of the plant (see Theory,
% e.g. phase change, etc.).

% In order to do ourself the pulling out of the deltas, we redefine the
% uncertainites in normalized terms, with a delta related to phi s.t.:
delta_phi_nom = 0;
delta_psi = ureal('delta_psi',delta_phi_nom,'range',[-1,1]);

% Indeed we need to use all deltas between -1 and 1, also as lftdata gives them to us.
% The weighting of these normalized deltas into the real uncertainty
% parameters is done via W (see later), as lftdata does building P.

% We can therefore express phi as: phi = delta*phi_max . 

% Coupling matrix C:
% C_mat = [ cos(phi) sin(phi); 
%          -sin(phi) cos(phi)];
% OBS.: these are uncertain parameters, not jw dependant because it is
% assumed that are constant in time, therefore in Laplace domain they
% remain unchanged.

% OBS.: nominal value of this coupling matrix is (for phi = 0) = identity
% matrix eye(2). So in the nominal case the controls entering the two 
% dynamics, u_lon and u_lat, are the ones actually produced by the regulators
% R_p and R_q, i.e.: u_lat = delta_lat and u_lat = delta_lat. Ok. 

% the uncertain variable psi is inside non-linear functions (-> lftdata won't work...TRY!), therefore when you
% would express them in terms of in-(-1)-to-1 range adimensional deltas as
% psi_bar*(1+r_psi*delta_psi), you would have your uncertainty in the
% diagram as blocks with inside cos(psi_bar*(1+r_psi*delta_psi)) and 
% sin(psi_bar*(1+r_psi*delta_psi)): you are unable to pull out the deltas
% to build the Delta matrix with only the deltas inside (as example in
% class with delta^2 etc.) (both manually and lftdata will be unable to
% (check this statement by running everything with the nonlinear
% uncertaintes - checked): a P-Delta form, which is used for exact uncertainty modelling,
%  cannot be achieved if the deltas are arguments of nonlinear functions.
%  Therefore the options we have are the following:
% MODEL MODEL UNCERTAINTESM - OPTIONS:
% OPTION 1 - STILL USING P-DELTA (meant to be exact) APPROACH
% 1) - two indipendent deltas: change of variables -> now the new uncertain
% variables are the 2 non-linear functions of delta: z = sin and c =cos, each in
% their possible range, given phi range and their multiplicative
% uncertainty models: z_nom*(1+r_z*delta_z) and z_nom*(1+r_z*delta_z). 
% ----> Downside: P-delta approach will loose its exact-uncertainty modelling
%                 property, as it is allowing for values of sin and cos
%                 that are meaningless, as they are indipendent
%                 uncertain/aleatory variables and do not need to satisfy
%                 sin^2 + cos^2 = 1, so you could even get c = 0.1 and
%                 sin=0.1. This method allows ALL possible relative values
%                 between c and s, therefore is the most conservative.
% 
% 2) - linearization: linearize the functions f_i(psi) in around the mean
% value of the uncertainty. Since psi aleatory variable is uniformely
% distributed in [-15°,15°], its mean value is 0°. Linearization in psi =
% 0° can be performed via McLaurin, giving: sin(psi) ~ psi , cos(psi) ~ 1.
% We stop to first order because the higher order terms would be non linear
% in psi, so non-linear in delta, re-introducing the problem: indeed we
% want to LINEARIZE. Observe that the new approximated f_i(psi) in both
% cases is overbounding the real f_i in our [-15°,15°] range: indeed you
% are always considering the cosine to be = 1 (highest), so this doesn't even have an
% associated delta anymore, and you are considering only the uncertainty in
% the sin: psi:
%   [ 1  psi;
%    -psi 1];
% Since the approximated functions are always higher in magnitude(!) than the real ones in
% the [-15°,15°] range, we are being conservative, but still considering
% all possible values of s and c (+ some non physical ones), indeed his 
% choice will generate combinations of s and c that are NEVER verifying the
% physical constraint s^2 + c^2 = 1 (unless in the case psi = 0, for which
% the uncertainty is not present). Notice that |f_approx_i| > |f_i| always 
% in our interval. In another scenario were the |f_i| were to be concave in
% psi_nom, the linearized |f_approx| (line tangent to f in psi_nom) would be
% below the real |f| in magnitude, and therefore by linearizing we would
% underestimate the uncertainty.
% 

% OPTION 2 - CIRCULAR OVERBOUND (meant to be approximate: conservative) approach
% You separate the 2x2 matrix in two blocks, one in cascade with G_lat(s)
%  and one in cascade with G_lon(s), and build the polar of this 2-blocks
%  cascade as uncertaint plants and see if a conservative circular
%  overbound (use ucover()...) is a good approximation or too restrictive. If the ucertainty
%  regions (4 of them, 1 for each of the elements of the 2x1 G_lat: G_p, G_p
%  and G_lon: G_q, G_theta) occupy a lot of the circular areas, then this
%  is a good approximation, otherwise we are being too conservative.
%  The additional regions I am including expanding the exact uncertainty
%  regions to circular ones are non-physical points (e.g. with sin^2 +
%  cos^2 not = 1), as for the case of the indipendent deltas approach.

% Choice of method: By using
% all three of them, if AT LEAST one of them satisfies the RS condition, then
% we satisfy it for a conservative, bigger uncertainty, so also the exact
% uncertain model would. If some of the three verify RS and some others do
% not, then the ones that do not are considering a too conservative
% approach. The first method is for sure the most conservative: the other 2
% take into account uncertainty couple values (c,s) that are a subset of
% this one, so if RS is verified for it, it will be for all others. At the
% same time, if RS is checked for another less conservative one, it won't
% be necessary check for uncertainty model 1. 

%%
% We will use...which? We try with the linearized one. Then, if it doesn't
% satisfy R.S., this means that it is either that our system is RS but our 
% uncertainty model is too conservative (considering too high
% uncertaintes), or that our system is actually not RS. If this happens, we
% try another of the two approaches, and then if the same happens we try
% with the third. If all do not verify RS, then we can only say that our
% system is not RS for a conservative uncertainty model, and we would need
% to find another way of modeling uncertainty as exact for the case of
% non-linear uncertainty relations. In this scenario only finding such 
% method (or another conservative (but less conservative) method) would 
% allows us to certainly say that our system is RS.

% OPTION 1.2 : linearize the functions of the uncertainties

% Coupling Block C:
% C_mat = [ 1   phi ; % Linearized Coupling matrix C
%          -phi  1 ];

% Pull out the deltas into a MULTIPLICATIVE FORM (possible given the
% linearity in delta):

% Weight
Wc_mat = [0 phi_max
          -phi_max, 0];
Wc_blk = ss(Wc_mat);  % block (uss: uncertain space state)
Wc_blk.u = {'delta_lon', 'delta_lat'};   % inputs to Wc
Wc_blk.y = {'z_C_lon', 'z_C_lat'};       % outputs of Wc

% Delta_c matrix
Delta_c_mat = eye(2)*delta_psi;       % Normalized in [-1,1]. Ok.
Delta_c = ss(Delta_c_mat);            
Delta_c.u = {'z_C_lon', 'z_C_lat'};   % inputs to Delta_c
Delta_c.y = {'w_C_lon', 'w_C_lat'};   % outputs of Delta_c

% Sum node: vectorial summing node -> as 2 sum nodes (one per component):
SumLong= sumblk('u_lon = delta_lon + w_C_lon');
SumLat = sumblk('u_lat = delta_lat + w_C_lat');

% System blocks
% The other blocks of each augmented system (lat and lon) are already inherited from 
% task1_lat,task2_lat,task3_lat (lateral dynamics closed loop system)
% and 
% task1_lon,task2_lon,task3_lon (longitudinal dynamics closed loop system).

% For each dynamics we consider the already decomposed dynamics of
% G_lat_unc, G_lon_unc into their G_lat, Delta_lat and G_lon, Delta_lon
% respectively (this has been done in task3_lat.m and task3_lon.m):
%           G_lat_unc -> G_lat, Delta_lat
%           G_lon_unc -> G_lon, Delta_lon
% Substitute nominal input delta to after-coupling input u (last elements of the lists):
G_lat.u = {'ud1_lat','ud2_lat','ud3_lat','ud4_lat','ud5_lat','u_lat'};
G_lon.u = {'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon','u_lon'};

% OBS.: 4x4 MIMO system

% Redefine inputs to the control limitation weight blocks W_q and to the plants:
W_q_lat.u = 'u_lat' ;
W_q_lon.u = 'u_lon';


% We do not include Delta_lat and Delta_lon, as we did in the task3 for the
% single dynamics.
M = connect(SumOuterLat,R_phi_tuned, R_p_tuned, G_lat, ... % ss
            Wc_blk, SumLat, SumLong, ... 
            SumOuterLon,R_theta_tuned,R_q_tuned, G_lon, ...
    {'ud1_lat','ud2_lat','ud3_lat','ud4_lat','ud5_lat', ... % inputs in M coming from Delta_lat
     'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon',...
     'w_C_lon', 'w_C_lat'}, ...
    {'yd1_lat','yd2_lat','yd3_lat','yd4_lat','yd5_lat', ... % outputs from M coming going into Delta_lon
     'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon',...
     'z_C_lat', 'z_C_lon'});   


%% 1) SIGMA_MAX of M
% M: 12x12 indeed we used 12 channels for Delta in connect(). 12>11 number
% of uncertain parametrs: ok. also in lecture (due to delta_psi acting on 2
% signals).
tf(M)
isstable(M) % check if M stable (necessary to check)
max(max(sigma(M))) % 74.5232


% Plot singular values of M(jw): sigma(M(jw))
wspan = logspace(-3, 4, 800);
[sv, wout] = sigma(M, wspan);      % rows = singular values, descending order
sv_dB = 20*log10(sv);

% colors for the plot
col_low = [0.80 0.82 0.85];    
col_bar = [0.42 0.56 0.68];      

figure; hold on; grid on; box on
set(gca,'XScale','log')

% lower sigma_i (rows 2..end)
for i = 2:size(sv_dB,1)
    semilogx(wout, sv_dB(i,:), '-', 'Color', col_low, 'LineWidth', 1.0);
end

% sigma_bar = first row (the highest): thick and highlighted in the figure
hbar = semilogx(wout, sv_dB(1,:), '-', 'Color', col_bar, 'LineWidth', 2.8);

yline(0, 'k:', 'LineWidth', 1, 'Label', '0 dB', 'LabelHorizontalAlignment','left');

% peak of sigma_bar, annotated in dB and linear
[pk_dB, ipk] = max(sv_dB(1,:));
w_pk = wout(ipk);
plot(w_pk, pk_dB, 'o', 'Color', col_bar, 'MarkerFaceColor', col_bar, 'MarkerSize', 7);
text(w_pk, pk_dB, sprintf('  %.1f dB  (%.1f lin)', pk_dB, 10^(pk_dB/20)), ...
     'Color', col_bar*0.7, 'FontSize', 9, 'VerticalAlignment','bottom');

% dummy handle for the legend of the lower sigmas
hlow = semilogx(nan, nan, '-', 'Color', col_low, 'LineWidth', 1.0);

xlabel('\omega  [rad/s]'); ylabel('Singular values  [dB]');
title('Singular values of M(j\omega)');
legend([hbar hlow], ...
    {'$\bar{\sigma}\,(M(j\omega))$', '$\sigma_i\,(M(j\omega)),\ i>1$'}, ...
    'Interpreter','latex', 'Location','southwest');
xlim([wspan(1) wspan(end)]);


%% 2) Scalings method
% To reduce conservativeness exploiting knowledge of Delta's diagonality.
% Choice of D(jw) s.t.: sigma_max((D^(-1)(jw))*M(jw)D(jw))< 1 for all w. 

% Methodology: M. K. H. Fan, A. L. Tits, J. C. Doyle (1991)
% "Robustness in the presence of mixed parametric uncertainty and unmodeled dynamics"
% IEEE Transactions on Automatic Control, vol. 36, no. 1, pp. 25–38.

% STEP 0 --- Problem setup -----------------------------------------------

% (A.) If we cycle on a vector of sampled w:
N = 800; % investigate 200 frequencies
wspan = logspace(-3, 4, N); % 10^-3 - 10^4
% Then, at each w=w_bar, the problem is numerical: (D is D(j*w_bar))

% D=diag(d1,...,d12), d_i > 0 s.t.sigma_max((D(jw))*M(jw)D^(-1)) < 1.
% where now sigma max is the heighest singular value.

% (B.) Normalize by setting d_12 = 1: indeed one dof (of 12) is redundant,
% since D M D^(-1) is invariant wrt global scaling of D (namely, using D or
% k*d is the same, so we fix the 'scale' od D by fixing one of the 12
% values of D already). We are left with 11 free variables: d_1, ..., d_11
% to use to minimize J = sigma_max(D M D^(-1))

% (C.) Operate a change of variables:
%  new free variables to minimize J: x_i = log(d_i) \in R   , i = 1,...,11
%  namely, the elements of D become d_i = e^(x_i). In this new variables
% x_i the problem is convex: only has one minimum, the global one ->
% assured convergence of the minimization! (Also said in the cited papers
% and in others)[Add rigorous proof?].

% STEP 1 -- Minimization for RS at a single frequency ---------------------
% For each w_k (k-th omega):
% -  use intial guess D0 = eye(12) (consistent with choice of d_12 = 1),
% meaning x0 = zeros(11,1) (only 11 variables: d_12 is fixed);
% -  evaluate J = sigma_max(M(j*w_k))
% -  minimize
% Until residual |J(t+1) - J(t)| < tol (tol: to be chosen)
% You have found minJ = min max sigma_max(DMD^(-1)) at this w_k. 
% -> if minJ < 1 ok, the syst. is robustly stable at THIS frequency. If not, then we can already discard the scaling method
% with this minimization approach and go to mu-analysis instead, as exists
% w s.t. sigma_max(D M D^(-1))>1, so we cannot say anything about R.S.
% x_i_opt -> get d_i opt = log(x_i_opt) at this frequency.

% STEP 2 -- Rational fit of D(jw) -----------------------------------------
%  for the check on whole frequency domain.
% We now have sigma_max() at only a bunch of frequencies (plot these
% points). We want to find the evolution of D(jw), therefore of d(jw) along
% w to have an approximated curve of sigma_max. Therefore 
% find d_i(s) rational, in s in general, such that have no zeros in RHP
% (because then they need to be inverted so RHP zeros would become RHP poles)
%, so minimum phase (otherwise D^(-1) unstable), and no RHP poles (otherwise D unstable)
% Since (each) d_i(s) = d_i(jw) will need ideally to pass through the
% d_i(w_k) values that are real in those wk, then there d_i(s) shall be
% (ideally) real, and therefore d_i(jw_k) = |d_i(jwk)|. So we basically
% find d(jw) by finding the TF s.t. its magnitude plot  |d_i(jwk)| fits the
% d_i(wk) dataset found, even though in those point its phase may not be =0,
% so d_i(jwk) will not be equale to the values found, but only their magnitude will. 
% -> find d_i(jw) , and CHECK on all 12 of them: no RHP poles (so stable) or RHP zeros (minimum phase).
% -> D(jw) is known.
% STEP 3 -- Verify RS condition -------------------------------------------
% - Compute and plot sigma_max (D(jw)*M(jw)*D(jw)^(-1))
% - Check that always <1 for all w. If not -> mu_analysis (see below).

d_opt = zeros(12,N); % allocate optimal d_i
M_tf = tf(M);
i = 0;
x0 = zeros(11,1);
d_max = 100;% bounds for optimization of d, that result in bounds on opt x:
x_max = log(d_max); 
x12 = 0; % fixes d_12 = 1.
for w = wspan
    i = i+1;
    M_w = evalfr(M_tf, 1j*w); % evaluate M(jw) in the current frequency w
    % Function to minimize (in x in R^11):
    J = @(x) max(svd( diag(exp([x;x12])) * M_w * diag(exp(-[x;x12])) )); 
    % max singular values of frequency response of M'
    lb = -x_max*ones(11,1); % limit values: 
    ub =  x_max*ones(11,1);
    x_opt = fmincon(J, x0, [],[],[],[],lb,ub,[], ...
                    optimoptions('fmincon','Display','off'));  
 % No need to use fmincon because constraint on d_i>0 sign descend
 % automatically from the change of variables we will impose to acquire d_i
 % from the optimum x_i.
  d_opt(:,i) = exp([x_opt;x12]);% since x real (because fmincon finds only
                                % real solutions), d_i is positive for sure.
    x0 = x_opt; % at next iteration start at previous minimum coordinate. 
end
%% Plot of descrete d_i : d11 to d12(=1)
figure
cmap = [
    0.65 0.35 0.35;  
    0.35 0.55 0.70;  
    0.40 0.68 0.55;   
    0.78 0.70 0.35;  
    0.70 0.50 0.35;  
    0.55 0.65 0.75;   
    0.60 0.40 0.65;   
    0.40 0.62 0.60;   
    0.75 0.45 0.45;   
    0.50 0.60 0.40;   
    0.65 0.55 0.70;  
    0.45 0.50 0.60;   
];

h = gobjects(12,1); % array of 12 graphic handles
for c = 1:12
    h(c) = loglog(wspan, d_opt(c,:)', ...
                  'Color', cmap(c,:), ...
                  'LineWidth', 1.8);
    hold on
end

grid on
xlabel('\omega  [rad/s]')
ylabel('d_i')
title('Optimal D-scaling coefficients  d_i(\omega)')

labels = arrayfun(@(k) sprintf('$d_{%d}$', k), 1:12, 'UniformOutput', false);
legend(h, labels, 'Interpreter', 'latex', 'Location', 'best', 'NumColumns', 2)

%% max sigma with descrete real scalings

sv_scaled = zeros(12, N);
for k = 1:N
    D_k = diag(d_opt(:,k));
    M_w = evalfr(M_tf, 1j*wspan(k));
    sv_scaled(:,k) = svd(D_k *M_w * eye(12)/D_k);
end
sv_scaled_dB = 20*log10(sv_scaled);
max(max(sv_scaled_dB));  % 12 sigmas that evolve in time
max(max(sv_scaled))

%% Plot of discretely scaled M:
%%Singular values of D^{-1}MD at the optimized discrete D

% Combined plot: unscaled vs scaled
figure; hold on; grid on; box on
set(gca,'XScale','log')

% --- Unscaled: lower sigmas ---
col_low_u  = [0.80 0.82 0.85];
col_bar_u  = [0.42 0.56 0.68];

for i = 2:size(sv_dB,1)
    semilogx(wout, sv_dB(i,:), '-', 'Color', col_low_u, 'LineWidth', 0.8);
end
hbar_u = semilogx(wout, sv_dB(1,:), '-', 'Color', col_bar_u, 'LineWidth', 2.5);

% --- Scaled: lower sigmas ---
col_low_s  = [0.85 0.80 0.78];
col_bar_s  = [0.72 0.38 0.35];

for i = 2:12
    semilogx(wspan, sv_scaled_dB(i,:), '-', 'Color', col_low_s, 'LineWidth', 0.8);
end
hbar_s = semilogx(wspan, sv_scaled_dB(1,:), '-', 'Color', col_bar_s, 'LineWidth', 2.5);

% --- 0 dB line ---
yline(0, 'k:', 'LineWidth', 1, 'Label', '0 dB', ...
      'LabelHorizontalAlignment','left');

% --- Annotate peaks ---
[pk_u, ipk_u] = max(sv_dB(1,:));
plot(wout(ipk_u), pk_u, 'o', 'Color', col_bar_u, ...
     'MarkerFaceColor', col_bar_u, 'MarkerSize', 7);
text(wout(ipk_u), pk_u, sprintf('  %.1f dB', pk_u), ...
     'Color', col_bar_u*0.7, 'FontSize', 9, 'VerticalAlignment','bottom');

[pk_s, ipk_s] = max(sv_scaled_dB(1,:));
plot(wspan(ipk_s), pk_s, 'o', 'Color', col_bar_s, ...
     'MarkerFaceColor', col_bar_s, 'MarkerSize', 7);
text(wspan(ipk_s), pk_s, sprintf('  %.1f dB', pk_s), ...
     'Color', col_bar_s*0.7, 'FontSize', 9, 'VerticalAlignment','bottom');

% --- Dummy handles for lower sigmas in legend ---
hlow_u = semilogx(nan, nan, '-', 'Color', col_low_u, 'LineWidth', 0.8);
hlow_s = semilogx(nan, nan, '-', 'Color', col_low_s, 'LineWidth', 0.8);

xlabel('\omega  [rad/s]'); 
ylabel('Singular values  [dB]');
title('Effect of D-scaling on singular values of M(j\omega)');
legend([hbar_u, hlow_u, hbar_s, hlow_s], ...
       {'$\bar{\sigma}(M)$',                '$\sigma_i(M),\ i>1$', ...
        '$\bar{\sigma}(DMD^{-1})$',          '$\sigma_i(DMD^{-1}),\ i>1$'}, ...
       'Interpreter','latex', 'Location','southwest');
xlim([wspan(1) wspan(end)]);

%%
% Fit d_i(s): use Robust Control Toolbox's fitmagfrd:
% frd/fitmagfrd - Fit frequency response magnitude data with minimum-phase state-space model using log-Chebyshev magnitude design
% This MATLAB function computes a stable, minimum-phase state-space model
%  of order N with a frequency response that approximates the
%   frequency-response model A.

n_order = 5; % choose order of fitting d_i(s)

d_tf = tf(zeros(1,12)); % tf vector

% figure
for i = 1:12
    frd_i = frd(d_opt(i,:), wspan);       
    d_tf(i) = fitmagfrd(frd_i, n_order);   % fit: stable + minimum phase: ok.
  %  bode(d_tf(i)); hold on; grid on;
end



% sigma scaled with complex scaling rational functions

sigma_opt = zeros(1,N);
for k = 1:N
    % Evaluate each d_i(s) in each s = jw_k 
    d_vals = zeros(12,1);
    for i = 1:12
        d_vals(i) = evalfr(d_tf(i), 1j*wspan(k));
    end
    D_k = diag(d_vals);          % D_k is now 12x12 complex
    
    M_w = evalfr(M_tf, 1j*wspan(k));
    sigma_opt(k) = max(svd(D_k * M_w * (eye(12)/D_k)));
end
max(sigma_opt) % 0.8721


%% Plot complex scalings

% --- compute sigma with fitted D(s) ---
sigma_fitted = zeros(1,N);
for k = 1:N
    d_vals = arrayfun(@(sys) evalfr(sys, 1j*wspan(k)), d_tf).';
    D_k = diag(d_vals);
    M_w = evalfr(M_tf, 1j*wspan(k));
    sigma_fitted(k) = max(svd(D_k * M_w * (eye(12)/D_k)));
end
fprintf('max sigma fitted: %.4f\n', max(sigma_fitted))

% --- colors (consistent with previous plots) ---
col_low_u = [0.80 0.82 0.85];   col_bar_u = [0.42 0.56 0.68];  % blue  (unscaled)
col_low_s = [0.85 0.80 0.78];   col_bar_s = [0.72 0.38 0.35];  % red   (discrete D)
col_bar_f = [0.35 0.62 0.40];                                  % green (fitted D(s))

fig_combined = figure; hold on; grid on; box on
set(gca,'XScale','log')

% --- 1) Unscaled ---
for i = 2:size(sv_dB,1)
    semilogx(wout, sv_dB(i,:), '-', 'Color', col_low_u, 'LineWidth', 0.8);
end
hbar_u = semilogx(wout, sv_dB(1,:), '-', 'Color', col_bar_u, 'LineWidth', 2.5);

% --- 2) Discrete D-scaling ---
for i = 2:12
    semilogx(wspan, sv_scaled_dB(i,:), '-', 'Color', col_low_s, 'LineWidth', 0.8);
end
hbar_s = semilogx(wspan, sv_scaled_dB(1,:), '-', 'Color', col_bar_s, 'LineWidth', 2.5);

% --- 3) Fitted D(s) ---
hbar_f = semilogx(wspan, 20*log10(sigma_fitted), '-', ...
                  'Color', col_bar_f, 'LineWidth', 2.5);

% --- 0 dB ---
yline(0, 'k:', 'LineWidth', 1, 'Label','0 dB', 'LabelHorizontalAlignment','left');

% --- peaks ---
[pk_u, ipk_u] = max(sv_dB(1,:));
plot(wout(ipk_u), pk_u, 'o', 'Color',col_bar_u, 'MarkerFaceColor',col_bar_u, 'MarkerSize',7);
text(wout(ipk_u), pk_u, sprintf('  %.1f dB',pk_u), ...
     'Color',col_bar_u*0.7, 'FontSize',9, 'VerticalAlignment','bottom');

[pk_s, ipk_s] = max(sv_scaled_dB(1,:));
plot(wspan(ipk_s), pk_s, 'o', 'Color',col_bar_s, 'MarkerFaceColor',col_bar_s, 'MarkerSize',7);
text(wspan(ipk_s), pk_s, sprintf('  %.1f dB',pk_s), ...
     'Color',col_bar_s*0.7, 'FontSize',9, 'VerticalAlignment','bottom');

[pk_f, ipk_f] = max(20*log10(sigma_fitted));
plot(wspan(ipk_f), pk_f, 'o', 'Color',col_bar_f, 'MarkerFaceColor',col_bar_f, 'MarkerSize',7);
text(wspan(ipk_f), pk_f, sprintf('  %.1f dB',pk_f), ...
     'Color',col_bar_f*0.7, 'FontSize',9, 'VerticalAlignment','bottom');

% --- dummy handles for minor sigmas ---
hlow_u = semilogx(nan,nan,'-','Color',col_low_u,'LineWidth',0.8);
hlow_s = semilogx(nan,nan,'-','Color',col_low_s,'LineWidth',0.8);

xlabel('\omega  [rad/s]');
ylabel('Singular values  [dB]');
title('RS verification — D-scaling methods comparison');
legend([hbar_u, hlow_u, hbar_s, hlow_s, hbar_f], ...
       {'$\bar{\sigma}(M)$',                                     '$\sigma_i(M),\ i>1$', ...
        '$\bar{\sigma}(D_{\rm opt}MD_{\rm opt}^{-1})$',          '$\sigma_i(D_{\rm opt}MD_{\rm opt}^{-1}),\ i>1$', ...
        '$\bar{\sigma}(D(j\omega)MD(j\omega)^{-1})$'}, ...
       'Interpreter','latex', 'Location','southwest');
xlim([wspan(1) wspan(end)]);

%%
figure
sigma_opt = zeros(1,N);
for k = 1:N
    w = wspan(k);
    D_k = diag(d_opt(:,k));
    M_w = evalfr(M_tf, 1j*w);
    sigma_opt(k) = max(svd(D_k *(M_w)*eye(12)/D_k));
end
loglog(wspan, sigma_opt); yline(1,'r--');
max(sigma_opt) % 0.8721



%% 3) Structural singular value: mu analysis
% State of the art is MATLAB Robust Control Toolbox's:
% mussv - Compute bounds on structured singular value (µ_Delta)
% upper bounds and lower bounds (current state of the art)

y_Delta = {'yd1_lat','yd2_lat','yd3_lat','yd4_lat','yd5_lat', ...
            'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon', ...
            'z_C_lon','z_C_lat'};                 % <-- lon, then lat
u_Delta  = {'ud1_lat','ud2_lat','ud3_lat','ud4_lat','ud5_lat', ...
            'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon', ...
            'w_C_lon','w_C_lat'};
Mmu = M(y_Delta, u_Delta);

% Structure of Delta block (consistent with signals defined here above):
%  - 10 scalars, REA, distinct  (5 Delta_lat + 5 Delta_lon)  ->  [-1 0]
%  (the minus sign indicates they are real deltas), 1 is th dimension of
%  each block, 0 means it is a diagonal block.
%  - delta_psi: scalar REPEATED twice (so NOT distinct uncertaintes: they are
% the SAME delta): so real 2x2, diagonal block, therefore its structure is 
% exoressed as: [-2 0]

Delta_struct = [ -ones(10,1), zeros(10,1) ;
                     -2     ,     0      ];

omega = logspace(-3, 4, 1000); % from 10^-3 to 10^4
[bounds,~] = mussv(frd(Mmu,omega), Delta_struct);

% Upper and lower bounds
ru = squeeze(bounds.ResponseData(1,1,:));   % upper bound
rl = squeeze(bounds.ResponseData(1,2,:));   % lower bound

% Peak of the upper bound
muPeak = max(ru);
fprintf('mu peak = %.4f  ->  RS %s\n', muPeak, ...
        string(muPeak<1));

%%
% plot mu_Delta bounds
figure
semilogx(bounds.Frequency, ru, bounds.Frequency, rl), grid on
legend('\mu upper bound','\mu lower bound')
xlabel('\omega [rad/s]'), ylabel('\mu'), title('RS - structured singular value')



%% PLOT TO SUM UP THE THREE METHODS

%% FINAL SUMMARY PLOT: sigma unscaled, D(s) scaling, mu bounds

col_bar_u  = [0.42 0.56 0.68];   %  blue   — unscaled sigma
col_bar_f  = [0.35 0.62 0.40];   %   — complex D(s) scaling
col_mu_up  = [0.60 0.40 0.65];   %  — mu upper bound
col_mu_lo  = [0.72 0.55 0.35];   %  — mu lower bound

figure; hold on; grid on; box on
set(gca, 'XScale', 'log')

%  1) sigma_max(M) unscaled
h1 = semilogx(wout, sv_dB(1,:), '-', 'Color', col_bar_u, 'LineWidth', 2.2);

%  2) sigma_max with fitted complex D(s)
h2 = semilogx(wspan, 20*log10(sigma_fitted), '-', 'Color', col_bar_f, 'LineWidth', 2.2);

%  3) mu upper bound 
h3 = semilogx(bounds.Frequency, 20*log10(ru), '--', 'Color', col_mu_up, 'LineWidth', 2.0);

%  4) mu lower bound 
h4 = semilogx(bounds.Frequency, 20*log10(rl), ':', 'Color', col_mu_lo,  'LineWidth', 2.0);

% 0 dB threshold 
yline(0, 'k:', 'LineWidth', 1.2, 'Label', '0 dB  (RS threshold)', ...
      'LabelHorizontalAlignment', 'left');

% annotate peaks 
[pk1, ip1] = max(sv_dB(1,:));
plot(wout(ip1), pk1, 'o', 'Color', col_bar_u, ...
     'MarkerFaceColor', col_bar_u, 'MarkerSize', 7);
text(wout(ip1), pk1, sprintf('  %.1f dB', pk1), ...
     'Color', col_bar_u*0.7, 'FontSize', 9, 'VerticalAlignment', 'bottom');

[pk2, ip2] = max(20*log10(sigma_fitted));
plot(wspan(ip2), pk2, 'o', 'Color', col_bar_f, ...
     'MarkerFaceColor', col_bar_f, 'MarkerSize', 7);
text(wspan(ip2), pk2, sprintf('  %.1f dB', pk2), ...
     'Color', col_bar_f*0.7, 'FontSize', 9, 'VerticalAlignment', 'bottom');

[pk3, ip3] = max(20*log10(ru));
plot(bounds.Frequency(ip3), pk3, 'o', 'Color', col_mu_up, ...
     'MarkerFaceColor', col_mu_up, 'MarkerSize', 7);
text(bounds.Frequency(ip3), pk3, sprintf('  %.1f dB', pk3), ...
     'Color', col_mu_up*0.7, 'FontSize', 9, 'VerticalAlignment', 'bottom');

xlabel('\omega  [rad/s]')
ylabel('Magnitude  [dB]')
title('RS verification — summary of methods')
legend([h1, h2, h3, h4], ...
       {'$\bar{\sigma}(M)$  — no scaling', ...
        '$\bar{\sigma}(D(j\omega)MD(j\omega)^{-1})$  — D-scaling', ...
        '$\bar{\mu}_\Delta$  — upper bound', ...
        '$\underline{\mu}_\Delta$  — lower bound'}, ...
       'Interpreter', 'latex', 'Location', 'southwest')
xlim([wspan(1) wspan(end)])