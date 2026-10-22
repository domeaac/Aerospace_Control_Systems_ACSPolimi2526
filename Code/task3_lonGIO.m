

%% Task 3

% clc
% clear
% close all

 % task1_lon_30;
 % task2_lon;
 % save('refWorkspace_t1t2_lon.mat')

 load('refWorkspace_t1t2_lon.mat')

%% Tuned Longitudinal Augmented system

% Build regulators blocks
% We extract them from the tuned longitudinal dynamics system

% Tuned Regulators:

% R_theta:
% 1x1: e_theta -> q_0
R_theta_tuned = getBlockValue(CL_lon, 'R_theta');  
R_theta_tuned.u = 'e_theta'; R_theta_tuned.y = 'q_0';

% R_q: 
% 2x1 complex matrix (fdt in MIMO framework!): {'q_0','q'} -> 'delta_lon'
R_q_tuned = getBlockValue(CL_lon, 'R_q');
R_q_tuned.u = {'q_0','q'}; R_q_tuned.y = 'delta_lon';

% Sum nodes
SumOuterLon = sumblk('e_theta = theta_0 - theta');

% Weights of the augmented system
% W_p_lon
% W_q_lon

% Uncertain Plant
G_lon_unc = sys_lon_unc;


%% Plant G_unc in P-Δ form/decomposition

% P stores the known plant G_lon model information
% Δ stores the uncertaintes of the plant G_lon

[G_lon,Delta_mat,BlockStruct,uncertaintyNorm] = lftdata(G_lon_unc);
% G_lon: ss (not uncertain); stores nominal info about G_lon.
% Delta_lon: umat

% build P block
G_lon.u = {'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon','delta_lon'};
G_lon.y = {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon','q','theta'};

% build Delta_lon block: 
% umat -> uss (uncertain system)
Delta_lon = ss(Delta_mat);   % ss, D=Delta_lon, A/B/C null
Delta_lon.u = {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon'};  % inputs of Delta_lon: y_Delta
Delta_lon.y = {'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon'};  % outputs of Delta_lon: u_Delta

% OBS.: at this point we could already compute M:
% Note: do not insert Delta_lon in the system to be connected 
M1 = connect(SumOuterLon,R_theta_tuned,R_q_tuned, G_lon, ...
    {'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon'},{'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon'}); % uss
% THIS IS CORRECT, indeed once you have the dynamics and the Delta_lon
% separated you can right away look for N11 = M which is the tf from the 
% u_delta to the y_delta only, in your system. Also by hand you could at
% this point already compute M, without passing through the equivalent P,
% K, and N and then extracting M from N (see this alternative procedure below).
%  THOUGH WE ARE NOT PASSING THROUGH THE ACTUAL P OF A P-DELTA
% FORM, BECAUSE A P-DELTA FORM HAS P FEEDBACKED BY AN EQUIVALENT K. And 
% this is instead what the assignment requests. So the previous P is actually only a G_nom. 
% The following code instead  computes M going into the theory steps of Skogestad.
% That will indeed find the same result, with the exceptions of numerical
% errors due to the difference in the two procedures.
figure
wspan = logspace(-5, 4, 500);
[sv1,wout1] = sigma(M1,wspan);
semilogx(wout1, 20*log10(sv1)', 'b-','LineWidth', 1.5);
maxsv1 = max(max(sv1)); % max: 0.9479 <1 (in linear scale): (ASSUMPTIONS...)
% -> FEEDBACK SYSTEM IS ROBUSTLY STABLE.

%% Channels dimensions

% Uncertainty channels: y_delta, u_delta
 dim_u_delta = length(Delta_lon.y);     % numero uscite di Delta_lon = righe di N11
 dim_y_delta = length(Delta_lon.u);     % numero ingressi di Delta_lon = colonne di N11

% Exogeneous inputs and outputs
dim_w = 1; % w = theta0
dim_z = 2; % z = {z1, z2}

% Equivalent controller input and output 
dim_v = 3; % input  to generalized K:  v = {theta0 theta q}
dim_u = 1; % output of generalized K:  u = delta_lon


%%  Find Generalized plant P

sys_P = connect(SumOuterLon,W_p_lon, W_q_lon,G_lon,... % only info about G_lon and W and not about R_i (ok: as in examples from Skogestad)
    {'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon','theta_0','delta_lon'}, ...                   % ingressi
    {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon','z_1','z_2','theta_0','theta','q'});            % uscite

P = tf(sys_P); % numerical. Indeed now all uncertainty is in Delta_lon. 10x1. 
% Indeed has outputs: 
% {y_delta, z, v} = {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon','z_1','z_2','theta_0','theta','q'}


%% Find Generalized controller K

sys_K = connect(SumOuterLon,R_theta_tuned,R_q_tuned,... % OBS: only R appears, no G_lon or W (ok, in accordance with examples from Skogestad)
    {'theta_0','theta','q'}, ...               % inputs : v vector
    'delta_lon');                         % outputs:   u vector (control variable)

K = tf(sys_K);  % 3x1: ok. Indeed: K: {v}->{u} is {theta_0,theta,q} -> {delta_lon}

%% Find N

n_v = length(K.u);   % ingressi di K   = uscite   di P, che K riceve (v)
n_u = length(K.y);   % uscite   di K   = ingressi di P, che K pilota (u)

sys_N = lft(sys_P, sys_K,n_u,n_v);   

N_dirty = tf(sys_N);

% Pulisci tutti gli elementi di sys_N
tol = 1e-10;

N_clean = pole_cleaning(N_dirty, tol);

N = N_clean; % 7x6: {u_delta, w} -> {y_delta, z} because u and v of K have 
% been lumped inside N by the lft (lower fractional transformation).
% Namely:
%{'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon', w} -> {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon', z}. ok

%% Extract M
% Estrai M = N11 (upper-left matrix block):
 M = N(1:dim_u_delta, 1:dim_y_delta); 
 % ok: it has 5 inputs and 5 outputs:
 %{'ud1_lon','ud2_lon','ud3_lon','ud4_lon','ud5_lon'} -> {'yd1_lon','yd2_lon','yd3_lon','yd4_lon','yd5_lon'}. ok,
 %it is the autonomous system we need for MIMO robust stability analysis.

M = pole_cleaning(M);
 % OBS.: even if I change output z = z1 z2 , e.g. in only z1, M is still
 % the same: ok.
%% singular values
isstable(M)
[sv, wout] = sigma(M, wspan);
hold on;
semilogx(wout, 20*log10(sv)', 'r--','LineWidth', 1.5); % plots all singular values
maxsv = max(max(sv)); % max: 0.9479 <1 (in linear scale): (ASSUMPTIONS...)
% -> FEEDBACK SYSTEM IS ROBUSTLY STABLE.
% - same results if not for numerical errors (sv ~ 10^(-400/20) = 10^(-20))
% - this proves that both the first method and the one following more
%   steps are analogous, as expected, and is allowing us to double check
%   our results.
% - OTHER METHODS CAN BE USED TO CHECK R.S..

%% Colors
c1  = [0.9 0.9 0.1];   
c2  = [0.2 0.5 0.9];
c1b = c1 * 0.60;
c2b = min(c2 * 0.65 + 0.10, 1);

%% Figure 1 – All SVs of M1 and M (dB, log w)
figure; hold on;

for i = 1:size(sv1,1)
    semilogx(wout1, 20*log10(sv1(i,:)), '-', 'Color', c1, 'LineWidth', 1.5);
end
for i = 1:size(sv,1)
    semilogx(wout, 20*log10(sv(i,:)), '--', 'Color', c2, 'LineWidth', 1.5);
end

h1 = semilogx(nan, nan, '-',  'Color', c1, 'LineWidth', 1.5);
h2 = semilogx(nan, nan, '--', 'Color', c2, 'LineWidth', 1.5);
legend([h1 h2], {'M_1 (direct)', 'M (via N_{11})'}, 'Location', 'southwest');

xlabel('Frequency (rad/s)'); ylabel('Singular Values (dB)');
title('Singular Values – M_1 vs M');
grid on; grid minor;
xlim([wspan(1) wspan(end)]);
set(gca, 'XScale', 'log');

%% Figure 2 – Top 2 SVs, peaks annotated for M only
sv1_top2 = sv1(1:2,:);
sv_top2  = sv(1:2,:);

figure; hold on;

h1a = semilogx(wout1, 20*log10(sv1_top2(1,:)), '-',  'Color', c1,  'LineWidth', 2);
h1b = semilogx(wout1, 20*log10(sv1_top2(2,:)), '-',  'Color', c1b, 'LineWidth', 1.5);
h2a = semilogx(wout,  20*log10(sv_top2(1,:)),  '--', 'Color', c2,  'LineWidth', 2);
h2b = semilogx(wout,  20*log10(sv_top2(2,:)),  '--', 'Color', c2b, 'LineWidth', 1.5);

yline(0, 'k:', 'LineWidth', 1, 'Label', '0 dB');

% Peaks for M1 (dot only, no dropline)
handles_M1 = [h1a, h1b];
for k = 1:2
    crv = sv1_top2(k,:);
    col = get(handles_M1(k), 'Color');
    [mv, ix] = max(crv);
    semilogx(wout1(ix), 20*log10(mv), 'o', ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 7);
end

% Peaks for M: dot + vertical dropline + w label on axis
handles_M = [h2a, h2b];
yl = ylim;
for k = 1:2
    crv     = sv_top2(k,:);
    col     = get(handles_M(k), 'Color');
    [mv, ix] = max(crv);
    wmax_k  = wout(ix);
    db_peak = 20*log10(mv);

    % dot at peak
    semilogx(wmax_k, db_peak, 's', ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 7);

    % dB label above peak
    text(wmax_k, db_peak, sprintf(' %.2f dB', db_peak), ...
        'Color', col, 'FontSize', 8, 'VerticalAlignment', 'bottom');

    % vertical dotted dropline
    semilogx([wmax_k wmax_k], [yl(1) db_peak], ':', ...
        'Color', col, 'LineWidth', 1.2);

    % w value on x-axis
    text(wmax_k, yl(1), sprintf('%.2e', wmax_k), ...
        'Color', col, 'FontSize', 8, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Rotation', 30);
end

legend([h1a h1b h2a h2b], ...
    {'\sigma_1(M_1)', '\sigma_2(M_1)', '\sigma_1(M)', '\sigma_2(M)'}, ...
    'Location', 'southwest');

xlabel('Frequency (rad/s)'); ylabel('Singular Values (dB)');
title('Top 2 Singular Values – M_1 vs M');
grid on; grid minor;
xlim([wspan(1) wspan(end)]);
set(gca, 'XScale', 'log');

 %% Functions

function sys_clean = pole_cleaning(sys, tol)
% POLE_CLEANING  Removes numerical noise from tf coefficients
%
% Inputs:
%   sys  : tf or ss system (SISO or MIMO)
%   tol  : absolute tolerance for zeroing small coefficients 
%          if not put as input, a default 1e-14 value is chosen
%
% Output:
%   sys_clean : cleaned tf system

if nargin < 2
    tol = 1e-14;   % default tolerance
end

sys_tf = tf(sys);
[ny, nu] = size(sys_tf);

N_cell = cell(ny, nu);
for i = 1:ny
    for j = 1:nu
        [num, den] = tfdata(sys_tf(i,j), 'v');
        num(abs(num) < tol) = 0;
        den(abs(den) < tol) = 0;
        N_cell{i,j} = tf(num, den);
    end
end

sys_clean = tf(zeros(ny, nu));
for i = 1:ny
    for j = 1:nu
        sys_clean(i,j) = N_cell{i,j};
    end
end

end

