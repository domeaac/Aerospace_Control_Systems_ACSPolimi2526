 clc ; clear ; close all
 PathsInit

%% Task 1: Nominal modelling and design

plot_nd = 1; % to plot nominal dynamics (nd) data

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - Longitudinal dynamics - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - (1.1) Plant modelling - - - - - - - - - - - - - - - - 

% extract ANT_X_data
Aircraft = ANT_X_data();
g = Aircraft.g;

Xu_nom = Aircraft.Stability.Xu;
Xq_nom = Aircraft.Stability.Xq;
Mq_nom = Aircraft.Stability.Mq;
Xd_nom = Aircraft.Control.Xd;
Md_nom = Aircraft.Control.Md;

% Nominal Longitudinal Dynamics State Space Model
A_lon_nom = [Xu_nom Xq_nom -g; 0 Mq_nom 0; 0 1 0];
B_lon_nom = [Xd_nom Md_nom 0]';
C_lon_nom = [0 1 0; 0 0 1];
D_lon_nom = [0 0]';

% Nominal Plant block
sys_lon_nom = ss(A_lon_nom,B_lon_nom,C_lon_nom,D_lon_nom, ...
    'StateName',  {'u (m/s)', 'q (rad/s)', 'theta (rad)'}, ... % state
    'InputName',  {'delta_lon'}, ... % input
    'OutputName', {'q', 'theta'}); % output

% - - - - - - - - - - (1.2) Plant analysis - - - - - - - - - - - - - - - - 

% ---- Single TFs G_q and G_theta, and their zeros and poles ----

sys_lon_nom_tf = tf(sys_lon_nom);

G_q_nom   = sys_lon_nom_tf(1); %  tf G_q : \delta_{lon}  -> q
G_theta_nom = sys_lon_nom_tf(2); %  tf G_theta: \delta_{lon} -> theta

% plot chunk 1 here

    z_q = tzero(sys_lon_nom(1)); % zeros of nom. G_q(s)
    p_q = pole(sys_lon_nom(1));  % poles of nom. G_q(s)

% plot chunk 2 here

    z_theta = tzero(sys_lon_nom(2)); % zeros of nom. G_theta(s)
    p_theta = pole(sys_lon_nom(2));  % poles of nom. G_theta(s)

% plot chunk 3 here

% Comment: we have stable or simply stable cancellations, therefore all
% hidden dynamics are not unstable. This results in:
% - G_q(s)  : no zeros; one unstable pole in 2.64.
% - G_theta(s): no zeros; one unstavle pole in 2.64 and one additional pole
%                       in the origin because theta is the integral of q.

% ----- Poles and zeros of plant G(s) -----

zeros_lon = tzero(sys_lon_nom);
poles_lon = pole(sys_lon_nom);

% plot chunk 4 here

% Comment: actually the zero and pole in -0.1068 are the ones that are
% cancelled within Gp and Gtheta, so after cancellation they do not appear in
% the total plant G. However, rather than computing tzero and poles of 
% minreal(sys_lon_nom) (system after cancellation), it is better so check
% that what has been cancelled is not critical, like in this case.

%  Pole-zero map

% plot chunk 5 here

% Comment: (plot is after cancellations) pole of G_theta is covering G_q's same pole,
% plus an additional one in s = 0.

% ------ Transfer functions ------

% Bode of Nominal Plants G_q(jw), G_theta(jw):

% plot chunk 6 here

% Comment: it is evident that G_theta(s) plot is the G_q(s)/s (type -1 due to additional pole).

% ------ Nominal plants step response -----

% plot chunk 7 here

% Comment: both dynamics are unstable, as expected by the RHP poles.
% Therefore a controller is needed.

%% - - - - - - - - (1.3) Feedback Attitude control - - - - - - - - - - - - -

%    - - - - - Fix controller structure/architecture: - - - -
%    - - - - - - - build parametric feedback system - - - - -

% Kinematics control block
R_theta = tunablePID('R_theta','P');
R_theta.u = 'e_theta'; % input name
R_theta.y = 'q_0';    % output name

% Dynamics control block
R_q = tunablePID2('R_q','PID');

% Fix the values that are not to be tuned:
R_q.c.Value = 0;    % fix value
R_q.c.Free = false; % it is not free

R_q.b.Value = 1;
R_q.b.Free = false;

R_q.Tf.Value = 0.02; % [s] time constant T_p of the filter pole. 
R_q.Tf.Free = false;

R_q.u = {'q_0','q'};      % inputs
R_q.y = {'delta_lon'}; % output (control variable)

% Outer summation node
SumOuter_lon = sumblk('e_theta = theta_0 - theta');

% Connect blocks
CL_tun = connect(SumOuter_lon,R_theta,R_q,sys_lon_nom,'theta_0',{'theta','q'});

%           - - - - Nominal design requirements - - - -

%%           ~ ~ ~ A) Nominal performance: W_p(jw) ~ ~ ~

% desired nominal dynamic performance
omega_n_lon = 10;       % [rad/s]
xi_lon = 0.9;

s = tf('s'); % define Laplace variable

F_2_theta = omega_n_lon^2/(s^2 + 2*s*xi_lon*omega_n_lon + omega_n_lon^2); % desired approx. F
S_2_theta = 1-F_2_theta;                                        % desired approx. S: constrained

% plot chunk 8 here

% Print step response info (T_set, S%, etc.)
F2_lon_stepresponse = stepinfo(F_2_theta);


% percentage overshoot and settling time

% Desired approx. F(jw)
% plot chunk 8b here

% Desired Nominal Sensitivity S_2(s) characteristics:

% 1) Steady state value (at w=0), = DC gain:
dc_lon = dcgain(S_2_theta);

% 2) Crossover frequency/bandwidth w_c:
[~, ~, ~, wc_lon] = margin(S_2_theta);   % Wcp = gain crossover

% 3) Resonance peak (damping)
[pk_lon, wpk_lon] = getPeakGain(S_2_theta);

% plot chunk 9 here

A_lon = 1e-3; % < A_d = -infty 
omega_b_lon = 6.71752; % > 0.9 rad/s, wrt w_cd = 7 rad/s
M_lon = 1.7165; % > M_d = 1.18
fprintf('omega_b_lon = %.4f rad/s\n', omega_b_lon);
fprintf('M_lon = %.4f\n', M_lon);

% Define performance weight W_p(s)

% METHOD 1 : manually (formulation above)
W_p_lon = (s/M_lon + omega_b_lon)/(s + A_lon*omega_b_lon);

% METHOD 2: via MATLAB function, with W of any order
W_p2inv_lon = makeweight(A_lon, omega_b_lon, M_lon, 0 , 1); 
% makeweight(A = DCgain, wc , M = HFgain, flag (0 = continuous time; 1 = discrete), order).

W_p2_lon = 1/W_p2inv_lon;

% comment to enforce method 1, uncomment to enforce method 2
% Plots show W_p to be a conservative requirement on S
W_p_lon = W_p2_lon;

% plot chunk 10 here

%%            ~ ~ ~ B) Control effort limitation: W_Q(jw) ~ ~ ~

alpha_lon = 5.1713; 
w_tau = 30;         

W_q_lon = alpha_lon*(s+w_tau*1e-3)/(s+w_tau);

% plot chunk 11 here

%         R.S: will be verified a posteriori
%         Design with structured Hinf synthesis

% Augmented plant blocks:
W_p_lon.u = 'e_theta';
W_p_lon.y = 'z_1';

W_q_lon.u = 'delta_lon';
W_q_lon.y = 'z_2';

% Augmented plant CL0 (Closed Loop iteration 0)
CL0_lon = connect(R_q, R_theta, sys_lon_nom, SumOuter_lon, W_p_lon, W_q_lon,... 
    'theta_0',{'z_1','z_2'},{'delta_lon','e_theta','theta','q_0','q'}); % input; output TO BE MINIMIZED; additional variables (u,e,y) (additional, not to be minimized, but useful later).

% options
opt_lon = hinfstructOptions('Display','final','RandomStart',20);
% solve design pb
[CL_lon,gamma_lon,info_lon] = hinfstruct(CL0_lon,opt_lon); % CL = solution CL with the solution/designed controllers; gamma=Hinf norm; info

%%                - - - - Analysis of the results - - - -

% a) Tuned/optimized Controllers parameters
if plot_nd
    showTunable(CL_lon)
end

% b) Nominal Complementary sensitivity
F_theta = getIOTransfer(CL_lon,'theta_0','theta'); % get Input Output tf that has input y0 = theta0 and y = output theta, therefore F(s) (nominal)

% pc 12 here
% Step response (of CL: F(s))
% pc 13 here
% pc 14 here

%% c) VERIFY PERFORMANCE 
S_lon = getIOTransfer(CL_lon,'theta_0','e_theta');

% pc 15 here
% d) VERIFY CONTROL LIMITATION in freq. domain

Q_lon = getIOTransfer(CL_lon,'theta_0','delta_lon');

% pc 16 here
% e) Nominal Loop TF L(s) : verify nominal wc.

L_lon = (1-S_lon)/S_lon;

% pc 17 here

% f) VERIFY CONTROL LIMITATION in time domain
t = linspace(0,6,10^4);
% doublet input
u = 0*(t<=1) + 10*(t>1 & t<= 3) - 10*(t>3 & t<= 5) + 0*(t>5); % step [°]
u = deg2rad(u); % [rad]

% pc 18 here

% (R.S. |M|_inf <1 will be verified in Task 2)

% 7) Verify Interpolation constraints (ICs)

% (R.P. |M|_inf + |Wp*S|_inf <1 will be verified in Task 4)

%% Plots

task1_lon_plots

fprintf('H-inf norm (gamma_lon) = %.4f  [target: < 1]\n', gamma_lon);


%% Summary of requirements vs obtained results 

% --- Desired 2nd order response ---
F2_info = stepinfo(F_2_theta);

% --- Obtained closed-loop response ---
F_theta_info = stepinfo(F_theta);

% --- Approx. to a second order via Giorgio's function ---
[wn_dom2, zeta_dom2, ~, ~] = compute_wc_xi(F_theta, linspace(0.001, 1000, 1000))

% --- Actual crossover frequency from S directly ---
[~, ~, ~, wc_cl] = margin(S_lon);

% --- Control effort: max delta_lon during doublet ---
delta_lon_sim = lsim(Q_lon, u, t);
max_delta_lon = max(abs(delta_lon_sim));

% PRINT SUMMARIES
fprintf('\n========================================\n')
fprintf('       REQUIREMENT A: NOMINAL PERFORMANCE\n')
fprintf('========================================\n')
pk_S2_lon = getPeakGain(S_2_theta);
pk_S_lon  = getPeakGain(S_lon);
fprintf('                  Required   Obtained - 2nd ord. approx    Obtained S\n')
fprintf('  wn   [rad/s] :  >= 10.00          %.4f                           -\n',  wn_dom2)
fprintf('  zeta         :  >= 0.90           %.4f                           -\n',  zeta_dom2)
fprintf('  wc   [rad/s] :  >= %.4f             -                         %.4f\n', wc_lon, wc_cl)
fprintf('  S peak  [-]  :  <= %.4f             -                         %.4f\n', pk_S2_lon, pk_S_lon)
if pk_S_lon <= pk_S2_lon && wc_cl >= wc_lon
    fprintf('  >> Requirement A: SATISFIED GLOBALLY (lower Speak, higher w_c than required)\n')
else
    fprintf('  >> Requirement A: NOT SATISFIED GLOBALLY\n')
    fprintf('  >> Requirement A: CHECK STEP RESPONSE INFO\n')
end

fprintf('\n========================================\n')
fprintf('       REQUIREMENT B: CONTROL EFFORT\n')
fprintf('========================================\n')
fprintf('                  Required     Obtained\n')
fprintf('  |delta_lon| :  <= 0.0500     %.4f\n', max_delta_lon)
if max_delta_lon <= 0.05
    fprintf('  >> Requirement B: SATISFIED\n')
else
    fprintf('  >> Requirement B: NOT SATISFIED\n')
end

fprintf('\n========================================\n')
fprintf('       STEP RESPONSE INFO\n')
fprintf('========================================\n')
fprintf('                  Desired 2nd order    Obtained\n')
fprintf('  Settling time [s] :   %.4f            %.4f\n', F2_info.SettlingTime,    F_theta_info.SettlingTime)
fprintf('  Rise time     [s] :   %.4f            %.4f\n', F2_info.RiseTime,        F_theta_info.RiseTime)
fprintf('  Peak value    [-] :   %.4f            %.4f\n', F2_info.Peak,            F_theta_info.Peak)
fprintf('  Overshoot     [%%] :   %.4f            %.4f\n', F2_info.Overshoot,       F_theta_info.Overshoot)
fprintf('========================================\n')

fprintf('\n  gamma (Hinf norm) : %.4f\n', gamma_lon)

%%

get_Ms_error = @(xi) M_lon - norm(tf([1, 2*xi*10, 0], [1, 2*xi*10, 10^2]), inf);
% 2. Solve for xi (using a reasonable guess below 0.707, since S only peaks if xi is low)
options = optimoptions('fsolve', 'Display', 'off');
new_required_xi_lon = fsolve(get_Ms_error, 0.4, options);

fprintf('Imposed S Peak (M_lon): %.4f\n', M_lon);
fprintf('Corresponding Requirement Closed-Loop xi: %.4f\n', new_required_xi_lon);

% Equivalent xi from the actual obtained S peak (same fsolve approach as sweep/func)
[pk_S_lon_actual, ~] = getPeakGain(S_lon);
get_Ms_error_actual  = @(xi) pk_S_lon_actual - norm(tf([1, 2*xi*10, 0], [1, 2*xi*10, 10^2]), inf);
xi_eq_from_pk_S      = fsolve(get_Ms_error_actual, 0.4, options);

fprintf('\nActual obtained S Peak:           %.4f\n', pk_S_lon_actual);
fprintf('Equivalent xi from actual S peak: %.4f\n', xi_eq_from_pk_S);


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
    magS = 20*log10(squeeze(magS));

end