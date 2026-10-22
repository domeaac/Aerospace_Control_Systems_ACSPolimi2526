 clc ;  clear ; close all
 PathsInit

%% Task 1: Nominal modelling and design

plot_nd = 1; % to plot nominal dynamics (nd) data

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - Lateral dynamics - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - (1.1) Plant modelling - - - - - - - - - - - - - - - - 

% extract ANT_X_data
Aircraft = ANT_X_data();
g = Aircraft.g;

Yv_nom = Aircraft.Stability.Yv;
Yp_nom = Aircraft.Stability.Yp;
Lp_nom = Aircraft.Stability.Lp;
Yd_nom = Aircraft.Control.Yd;
Ld_nom = Aircraft.Control.Ld;


% Nominal Lateral Dynamics State Space Model
A_lat_nom = [Yv_nom Yp_nom g; 0 Lp_nom 0; 0 1 0];
B_lat_nom = [Yd_nom Ld_nom 0]';
C_lat_nom = [0 1 0; 0 0 1];
D_lat_nom = [0 0]';

% Nominal Plant block
sys_lat_nom = ss(A_lat_nom,B_lat_nom,C_lat_nom,D_lat_nom, ...
    'StateName',  {'v (m/s)', 'p (rad/s)', 'phi (rad)'}, ... % state
    'InputName',  {'delta_lat'}, ... % input
    'OutputName', {'p', 'phi'}); % output

% - - - - - - - - - - (1.2) Plant analysis - - - - - - - - - - - - - - - - 

% ---- Single TFs G_p and G_phi, and their zeros and poles ----
sys_lat_nom_tf = tf(sys_lat_nom);

G_p_nom   = sys_lat_nom_tf(1); %  tf G_p : \delta_{lat}  -> p
G_phi_nom = sys_lat_nom_tf(2); %  tf G_phi: \delta_{lat} -> phi

% plot chunk 1 here

    z_p = tzero(sys_lat_nom(1)); % zeros of nom. G_p(s)
    p_p = pole(sys_lat_nom(1));  % poles of nom. G_p(s)

% plot chunk 2 here

    z_phi = tzero(sys_lat_nom(2)); % zeros of nom. G_phi(s)
    p_phi = pole(sys_lat_nom(2));  % poles of nom. G_phi(s)

% plot chunk 3 here

% Comment: we have stable or simply stable cancellations, therefore all
% hidden dynamics are not unstable. This results in:
% - G_p(s)  : no zeros; one unstable pole in 2.64.
% - G_phi(s): no zeros; one unstable pole in 2.64 and one additional pole
%                       in the origin because phi is the integral of p.

% ----- Poles and zeros of plant G(s) -----

zeros_lat = tzero(sys_lat_nom);
poles_lat = pole(sys_lat_nom);

% plot chunk 4 here

% Comment: actually the zero and pole in -0.1068 are the ones that are
% cancelled within Gp and Gphi, so after cancellation they do not appear in
% the total plant G. However, rather than computing tzero and poles of 
% minreal(sys_lat_nom) (system after cancellation), it is better so check
% that what has been cancelled is not critical, like in this case.

%  Pole-zero map

% plot chunk 5 here
% Comment: (plot is after cancellations) pole of G_phi is covering G_p's same pole,
% plus an additional one in s = 0.

% ------ Transfer functions ------

% Bode of Nominal Plants G_p(jw), G_phi(jw):

% plot chunk 6 here

% Comment: it is evident that G_phi(s) plot is the G_p(s)/s (type -1 due to additional pole).

% ------ Nominal plants step response -----

% plot chunk 7 here

% Comment: both dynamics are unstable, as expected by the RHP poles.
% Therefore a controller is needed.

%% - - - - - - - - (1.3) Feedback Attitude control - - - - - - - - - - - - -

%    - - - - - Fix controller structure/architecture: - - - -
%    - - - - - - - build parametric feedback system - - - - -

% Kinematics control block
R_phi = tunablePID('R_phi','P');
R_phi.u = 'e_phi'; % input name
R_phi.y = 'p_0';    % output name

% Dynamics control block
R_p = tunablePID2('R_p','PID');

% Fix the values that are not to be tuned:
R_p.c.Value = 0;    % fix value
R_p.c.Free = false; % it is not free

R_p.b.Value = 1;
R_p.b.Free = false;

R_p.Tf.Value = 0.02; % [s] time constant T_p of the filter pole. % ERA 0.02
R_p.Tf.Free = false;

R_p.u = {'p_0','p'};      % inputs
R_p.y = {'delta_lat'};    % output (control variable)

% Outer summation node
SumOuter_lat = sumblk('e_phi = phi_0 - phi');

% Connect blocks
CL_tun = connect(SumOuter_lat,R_phi,R_p,sys_lat_nom,'phi_0',{'phi','p'});

%           - - - - Nominal design requirements - - - -

%%           ~ ~ ~ A) Nominal performance: W_p(jw) ~ ~ ~

% desired nominal dynamic performance
omega_n_lat = 10;       % [rad/s] % 5 works
xi_lat = 0.9;

s = tf('s'); % define Laplace variable

F_2_phi = omega_n_lat^2/(s^2 + 2*s*xi_lat*omega_n_lat + omega_n_lat^2); % desired approx. F
S_2_phi = 1-F_2_phi;                                        % desired approx. S: constrained

% plot chunk 8 here

% Print step response info (T_set, S%, etc.)
F2_lat_stepresponse = stepinfo(F_2_phi);
% percentage overshoot and settling time

% Desired approx. F(jw)

% plot chunk 8b here

% Desired Nominal Sensitivity S_2(s) characteristics:

% 1) Steady state value (at w=0), = DC gain:
dc_lat = dcgain(S_2_phi);

% 2) crossover frequency wc_lat :
[~, ~, ~, wc_lat] = margin(S_2_phi);   % Wcp = gain crossover

% 3) peak
[pk_lat, wpk_lat] = getPeakGain(S_2_phi);

% plot chunk 9 here

A_lat = 1e-3; % < A_d = -infty 
%omega_b_lat = wc_lat;

omega_b_lat = 7.070107; %  9 rad/s, wrt w_cd = 7 rad/s
M_lat = 1.65731; %1.45*pk_lat; % > M_d = 1.18
fprintf('omega_b_lat = %.4f rad/s\n', omega_b_lat);
fprintf('M_lat = %.4f\n', M_lat);

% Define performance weight W_p(s)

% METHOD 1 : manually (formulation above)
W_p_lat = (s/M_lat + omega_b_lat)/(s + A_lat*omega_b_lat);
% METHOD 2: via MATLAB function, with W of any order 
W_p2inv_lat = makeweight(A_lat, omega_b_lat, M_lat, 0 , 1); 
% makeweight(A = DCgain, wc , M = HFgain, flag (0 = continuous time; 1 = discrete), order).

W_p2_lat = 1/W_p2inv_lat;

% comment to enforce method 1, uncomment to enforce method 2
% Plots show W_p to be a conservative requirement on S
W_p_lat = W_p2_lat;

% plot chunk 10 here

%%            ~ ~ ~ B) Control effort limitation: W_Q(jw) ~ ~ ~

alpha_lat = 4.65421; % bounded control action
w_tau = 30;                     % actuator speed % ERA 10!!!
                              % 50 rad/s implies a little less than 10 Hz
                              % it should be achievable by quadrotors actuator dynamics 

W_q_lat = alpha_lat*(s+w_tau*1e-3)/(s+w_tau);

% plot chunk 11 here

%         R.S: will be verified a posteriori 
%         Design with structured Hinf synthesis

% Augmented plant blocks:
W_p_lat.u = 'e_phi';
W_p_lat.y = 'z_1';
W_q_lat.u = 'delta_lat';
W_q_lat.y = 'z_2';

% Augmented plant CL0 (Closed Loop iteration 0)
CL0_lat = connect(R_p, R_phi, sys_lat_nom, SumOuter_lat, W_p_lat, W_q_lat,... 
    'phi_0',{'z_1','z_2'},{'delta_lat','e_phi','phi','p_0', 'p'}); % input; output TO BE MINIMIZED; additional variables (u,e,y) (additional, not to be minimized, but useful later).

% options
opt_lat = hinfstructOptions('Display','final','RandomStart',20);
% solve design pb
[CL_lat,gamma_lat,info_lat] = hinfstruct(CL0_lat,opt_lat); % CL = solution CL with the solution/designed controllers; gamma=Hinf norm; info

%%                - - - - Analysis of the results - - - -

% a) Tuned/optimized Controllers parameters
if plot_nd
    showTunable(CL_lat)
end

% b) Nominal Complementary sensitivity
F_phi = getIOTransfer(CL_lat,'phi_0','phi'); % get Input Output tf that has input y0 = phi0 and y = output phi, therefore F(s) (nominal)

% pc 12 here
% Step response (of CL: F(s))
% pc 13 here
% pc 14 here

%% c) VERIFY PERFORMANCE 
S_lat = getIOTransfer(CL_lat,'phi_0','e_phi');

% pc 15 here
% d) VERIFY CONTROL LIMITATION in freq. domain

Q_lat = getIOTransfer(CL_lat,'phi_0','delta_lat');

% pc 16 here
% e) Nominal Loop TF L(s) : verify nominal wc.

L_lat = (1-S_lat)/S_lat;

% pc 17 here

% f) VERIFY CONTROL LIMITATION in time domain
t = linspace(0,6,10^4);
% doublet input:
u = 0*(t<=1) + 10*(t>1 & t<= 3) - 10*(t>3 & t<= 5) + 0*(t>5); % step [°]
u = deg2rad(u); % [rad]

% pc 18 here

% (R.S. |M|_inf <1 will be verified in Task 2)

% 7) Verify Interpolation constraints (ICs)
% (R.P. |M|_inf + |Wp*S|_inf <1 will be verified in Task 4)

%% Plots

task1_lat_plots

fprintf('H-inf norm (gamma_lat) = %.4f  [target: < 1]\n', gamma_lat);


%% Summary of requirements vs obtained results

% --- Desired 2nd order response ---
F2_info = stepinfo(F_2_phi);

% --- Obtained closed-loop response ---
F_phi_info = stepinfo(F_phi);

% --- Approx. to a second order via Giorgio's function ---
[wn_dom2, zeta_dom2, ~, ~] = compute_wc_xi(F_phi, linspace(0.001, 1000, 1000))

% --- Actual crossover frequency from F directly ---
[~, ~, ~, wc_cl] = margin(S_lat)

% --- Control effort: max delta_lat during doublet ---
delta_lat_sim = lsim(Q_lat, u, t);
max_delta_lat = max(abs(delta_lat_sim));

% PRINT SUMMARIES
fprintf('\n========================================\n')
fprintf('       REQUIREMENT A: NOMINAL PERFORMANCE\n')
fprintf('========================================\n')
pk_S2_lat = getPeakGain(S_2_phi);
pk_S_lat  = getPeakGain(S_lat);
fprintf('                  Required   Obtained - 2nd ord. approx    Obtained S\n')
fprintf('  wn   [rad/s] :  >= 10.00          %.4f                           -\n',  wn_dom2)
fprintf('  zeta         :  >= 0.90           %.4f                           -\n',  zeta_dom2)
fprintf('  wc   [rad/s] :  >= %.4f             -                         %.4f\n', wc_lat, wc_cl)
fprintf('  S peak  [-]  :  <= %.4f             -                         %.4f\n', pk_S2_lat, pk_S_lat)
if pk_S_lat <= pk_S2_lat && wc_cl >= wc_lat
    fprintf('  >> Requirement A: SATISFIED GLOBALLY (lower Speak, higher w_c than required)\n')
else
    fprintf('  >> Requirement A: NOT SATISFIED GLOBALLY\n')
    fprintf('  >> Requirement A: CHECK STEP RESPONSE INFO\n')
end

fprintf('\n========================================\n')
fprintf('       REQUIREMENT B: CONTROL EFFORT\n')
fprintf('========================================\n')
fprintf('                  Required     Obtained\n')
fprintf('  |delta_lat| :  <= 0.0500     %.4f\n', max_delta_lat)
if max_delta_lat <= 0.05
    fprintf('  >> Requirement B: SATISFIED\n')
else
    fprintf('  >> Requirement B: NOT SATISFIED\n')
end

fprintf('\n========================================\n')
fprintf('       STEP RESPONSE INFO\n')
fprintf('========================================\n')
fprintf('                  Desired 2nd order    Obtained\n')
fprintf('  Settling time [s] :   %.4f            %.4f\n', F2_info.SettlingTime,  F_phi_info.SettlingTime)
fprintf('  Rise time     [s] :   %.4f            %.4f\n', F2_info.RiseTime,      F_phi_info.RiseTime)
fprintf('  Peak value    [-] :   %.4f            %.4f\n', F2_info.Peak,          F_phi_info.Peak)
fprintf('  Overshoot     [%%] :   %.4f            %.4f\n', F2_info.Overshoot,     F_phi_info.Overshoot)
fprintf('========================================\n')

fprintf('\n  gamma (Hinf norm) : %.4f\n', gamma_lat)

%% 

% compute and explain M requirement in terms of an equivalent xi 
get_Ms_error = @(xi) M_lat - norm(tf([1, 2*xi*10, 0], [1, 2*xi*10, 10^2]), inf);
% 2. Solve for xi (using a reasonable guess below 0.707, since S only peaks if xi is low)
options = optimoptions('fsolve', 'Display', 'off');
new_required_xi_lat = fsolve(get_Ms_error, 0.4, options);

fprintf('Imposed S Peak (M_lat): %.4f\n', M_lat);
fprintf('Corresponding Requirement Closed-Loop xi: %.4f\n', new_required_xi_lat);

% Equivalent xi from the actual obtained S peak (same fsolve approach as sweep/func)
[pk_S_lat_actual, ~] = getPeakGain(S_lat);
get_Ms_error_actual  = @(xi) pk_S_lat_actual - norm(tf([1, 2*xi*10, 0], [1, 2*xi*10, 10^2]), inf);
xi_eq_from_pk_S      = fsolve(get_Ms_error_actual, 0.4, options);

fprintf('\nActual obtained S Peak:           %.4f\n', pk_S_lat_actual);
fprintf('Equivalent xi from actual S peak: %.4f\n', xi_eq_from_pk_S);



% Giorgios funct: Define a function that given F, computes wc and xi:

function [wc, xi, magF, magS] = compute_wc_xi(F, omega)

% COMPUTE_WC_XI  Computes crossover frequency, damping ratio and phase margin
%                from complementary sensitivity function F (SISO only)
%
% Inputs:
%   F     : complementary sensitivity function (tf or ss, SISO)
%
% Outputs:
%   wc    : crossover frequency [rad/s] s.t. |L|=1
%   xi    : damping ratio (2nd order approximation)
%   phi_m : phase margin [deg]

    % --- Recover loop gain from F ---
    L = F / (1 - F); % from L def... is it FORMALLY correct to do it like this in MIMO considering a single F of whole F matrix? What exactly is this L?

    % --- Bode of L ---
    [mag, phase, wout] = bode(L, omega);
    mag   = 20*log10(squeeze(mag));   % [dB]
    phase = squeeze(phase);           % [deg] squeeze -> size: [1 x 1 x N] -> [N x 1]
    wout  = squeeze(wout);                  % squeeze -> size: [1 x 1 x N] -> [N x 1]

    % --- A.1 Crossover frequency: |L|_dB = 0 dB ---
    wc = interp1(mag, wout, 0);

    % --- A.2 Crossover phase and phase margin ---
    % Since F is not necessary of order 2, the definition of xi is
    % through its 2nd order representation:
    phi_c = interp1(wout, phase, wc); % crossover phase (reverse axis: find/interpolate phase_i at which wout_i = wc_i)
    phi_m = 180 - abs(phi_c);         % phase margin [deg]
    xi = sin(deg2rad(phi_m) / 2);     % exact formula for 2nd order approximation of F
    % alternatively: xi ≈ phi_m / 100 
     
    % F bode mag plot
    [magF, ~, ~] = bode(F, omega);
    magF = 20*log10(magF);

    % S bode mag plot
    S = 1 - F;
    [magS, ~, ~] = bode(S, omega);
    magS = 20*log10(magS);

end