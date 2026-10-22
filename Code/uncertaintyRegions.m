clc
clear
close all

%% Uncertainty regions

% We do this to see if the model uncertainty can be approximated using
% circular overbounding without introducing too much conservatorism (we
% would verify RS for other plants wrt the expected ones too) or not, as an
% alternative to the exact uncertain modelling (P-Δ): we have 11 uncertain
% parameters, so lftdata (symbolic) still works ok, but if the uncertainty
% regions are almost circular, that approach can be still a valuable
% alternative. If not, we have to rely necessarily on the exact one (P-Δ).

% COMMENT on results of the plots: since approximating the uncertainty
% regions with circular overbounds would introduce to much conservative
% robustness, we opt fot the exact P-Δ analysis, also because we expect it
% to behave well seing that we are only considering 5 parameters to be
% uncertain for each dynamics (namely: lateral and longitudinal).

% - - - - - - - Lateral Dynamics Plant uncertainty set Π - - - - - - - - -

% SIMO LTI system 

Aircraft = ANT_X_data();
g = Aircraft.g;

Yv_nom = Aircraft.Stability.Yv;
Yp_nom = Aircraft.Stability.Yp;
Lp_nom = Aircraft.Stability.Lp;
Yd_nom = Aircraft.Control.Yd;
Ld_nom = Aircraft.Control.Ld;

% uncertainties: +- 3 sigma as deterministic bounds
% (Uncertainties are given in terms of standard deviations)

unc_Yv = 3 * Aircraft.Stability.Yv_sigma_pct;   
unc_Yp = 3 * Aircraft.Stability.Yp_sigma_pct;   
unc_Lp = 3 * Aircraft.Stability.Lp_sigma_pct;   
unc_Yd = 3 * Aircraft.Control.Yd_sigma_pct;    
unc_Ld = 3 * Aircraft.Control.Ld_sigma_pct;     

fprintf('Uncertainty bounds \n')
fprintf('Yv: ±%.2f %%\n', unc_Yv)
fprintf('Yp: ±%.2f %%\n', unc_Yp)
fprintf('Lp: ±%.2f %%\n', unc_Lp)
fprintf('Yd: ±%.2f %%\n', unc_Yd)
fprintf('Ld: ±%.2f %%\n', unc_Ld)

% uncertain parameters
Yv_real = ureal('Yv_real', Yv_nom, 'Perc', unc_Yv);
Yp_real = ureal('Yp_real', Yp_nom, 'Perc', unc_Yp);
Lp_real = ureal('Lp_real', Lp_nom, 'Perc', unc_Lp);
Yd_real = ureal('Yd_real', Yd_nom, 'Perc', unc_Yd);
Ld_real = ureal('Ld_real', Ld_nom, 'Perc', unc_Ld);

% Uncertain State-Space Model
% Uncertain Lateral Dynamics State Space Model
A_lat_unc = [Yv_real Yp_real g; 0 Lp_real 0; 0 1 0];
B_lat_unc = [Yd_real Ld_real 0]';
C_lat_unc = [0 1 0; 0 0 1];
D_lat_unc = [0 0]';


% Uncertain system
G_lat_unc = ss(A_lat_unc, B_lat_unc, C_lat_unc, D_lat_unc); 
G_lat_unc.u = {'delta_lat'};
G_lat_unc.y = {'p', 'phi'};

% Nominal system
G_lat_nom = G_lat_unc.NominalValue;
G_lat_nom.u = {'delta_lat'};
G_lat_nom.y = {'p', 'phi'};


% Compute realizations

% p, phi tfs:

N = 100; % population

% G_p
G_lat_p_unc_ARRAY = usample(G_lat_unc(1),N);
G_lat_p_nom = G_lat_nom(1);
zpk(G_lat_p_nom)

% G_phi
G_lat_phi_unc_ARRAY = usample(G_lat_unc(2),N);
G_lat_phi_nom = G_lat_nom(2);
zpk(G_lat_phi_nom)

 
omega_full    = logspace(-2, 3, 500);
omega_regions = logspace(-1, 2, 8);

% G_p
uncertaintyRegionsPlot(G_lat_p_nom, G_lat_p_unc_ARRAY, N, ...
    omega_full, omega_regions, ...
    'Uncertainty Regions: G_{lat,p}(j\omega) — \delta_{lat} \rightarrow p')

% G_phi
uncertaintyRegionsPlot(G_lat_phi_nom, G_lat_phi_unc_ARRAY, N, ...
    omega_full, omega_regions, ...
    'Uncertainty Regions: G_{lat,\phi}(j\omega) — \delta_{lat} \rightarrow \phi')



%% Plots

omega = logspace(-1, 3, 500); % 500 points logaritmically spaced in frequency interval (10^(-2) - 10^3: [0.01, ...,1000]

% --------------------------- Gp(s)---------------------------------------
% G_p — Bode
plot_bode_unc(G_lat_p_nom, G_lat_p_unc_ARRAY, N, omega, ...
    'Bode: G_{lat,p}(j\omega) — \delta_{lat} \rightarrow p', ...
    '\omega [rad/s]', 'Magnitude [dB]', 'Nominal', 'Uncertain samples')

% G_p — Polar
plot_polar_uncertain(G_lat_p_nom, G_lat_p_unc_ARRAY, N, omega, ...
    'Polar: G_{lat,p}(j\omega) — \delta_{lat} \rightarrow p', ...
    'Nominal', 'Uncertain samples')


% --------------------------- Gphi(s) -------------------------------------

% G_phi — Bode
plot_bode_unc(G_lat_phi_nom, G_lat_phi_unc_ARRAY, N, omega, ...
    'Bode: G_{lat,\phi}(j\omega) — \delta_{lat} \rightarrow \phi', ...
    '\omega [rad/s]', 'Magnitude [dB]', 'Nominal', 'Uncertain samples')

% G_phi — Polar
plot_polar_uncertain(G_lat_phi_nom, G_lat_phi_unc_ARRAY, N, omega, ...
    'Polar: G_{lat,\phi}(j\omega) — \delta_{lat} \rightarrow \phi', ...
    'Nominal', 'Uncertain samples')



% - - - - - - - Longitudinal Dynamics Plant uncertainty set Π - - - - - - -


% SIMO LTI system

Xu_nom = Aircraft.Stability.Xu;
Xq_nom = Aircraft.Stability.Xq;
Mq_nom = Aircraft.Stability.Mq;
Xd_nom = Aircraft.Control.Xd;
Md_nom = Aircraft.Control.Md;

% uncertainties: +- 3 sigma as deterministic bounds
unc_Xu = 3 * Aircraft.Stability.Xu_sigma_pct;
unc_Xq = 3 * Aircraft.Stability.Xq_sigma_pct;
unc_Mq = 3 * Aircraft.Stability.Mq_sigma_pct;
unc_Xd = 3 * Aircraft.Control.Xd_sigma_pct;
unc_Md = 3 * Aircraft.Control.Md_sigma_pct;

fprintf('Longitudinal uncertainty bounds \n')
fprintf('Xu: ±%.2f %%\n', unc_Xu)
fprintf('Xq: ±%.2f %%\n', unc_Xq)
fprintf('Mq: ±%.2f %%\n', unc_Mq)
fprintf('Xd: ±%.2f %%\n', unc_Xd)
fprintf('Md: ±%.2f %%\n', unc_Md)

% uncertain parameters
Xu_real = ureal('Xu_real', Xu_nom, 'Perc', unc_Xu);
Xq_real = ureal('Xq_real', Xq_nom, 'Perc', unc_Xq);
Mq_real = ureal('Mq_real', Mq_nom, 'Perc', unc_Mq);
Xd_real = ureal('Xd_real', Xd_nom, 'Perc', unc_Xd);
Md_real = ureal('Md_real', Md_nom, 'Perc', unc_Md);

% Uncertain State-Space Model
% Uncertain Longitudinal Dynamics State Space Model
A_lon_unc = [Xu_real Xq_real -g; 0 Mq_real 0; 0 1 0];
B_lon_unc = [Xd_real Md_real 0]';
C_lon_unc = [0 1 0; 0 0 1];
D_lon_unc = [0 0]';

% Uncertain system
G_lon_unc = ss(A_lon_unc, B_lon_unc, C_lon_unc, D_lon_unc);
G_lon_unc.u = {'delta_lon'};
G_lon_unc.y = {'q', 'theta'};

% Nominal system
G_lon_nom = G_lon_unc.NominalValue;
G_lon_nom.u = {'delta_lon'};
G_lon_nom.y = {'q', 'theta'};

% Compute realizations
N = 500;

% G_q
G_lon_q_unc_ARRAY = usample(G_lon_unc(1), N);
G_lon_q_nom       = G_lon_nom(1);
zpk(G_lon_q_nom)

% G_theta
G_lon_theta_unc_ARRAY = usample(G_lon_unc(2), N);
G_lon_theta_nom       = G_lon_nom(2);
zpk(G_lon_theta_nom)

% Uncertainty regions
omega_full    = logspace(-2, 3, 500);
omega_regions = logspace(-1, 2, 8);

% G_q
uncertaintyRegionsPlot(G_lon_q_nom, G_lon_q_unc_ARRAY, N, ...
    omega_full, omega_regions, ...
    'Uncertainty Regions: G_{lon,q}(j\omega) — \delta_{lon} \rightarrow q')

% G_theta
uncertaintyRegionsPlot(G_lon_theta_nom, G_lon_theta_unc_ARRAY, N, ...
    omega_full, omega_regions, ...
    'Uncertainty Regions: G_{lon,\theta}(j\omega) — \delta_{lon} \rightarrow \theta')

% Plots
omega = logspace(-1, 3, 500);

% --------------------------- Gq(s) --------------------------------------
% G_q — Bode
plot_bode_unc(G_lon_q_nom, G_lon_q_unc_ARRAY, N, omega, ...
    'Bode: G_{lon,q}(j\omega) — \delta_{lon} \rightarrow q', ...
    '\omega [rad/s]', 'Magnitude [dB]', 'Nominal', 'Uncertain samples')

% G_q — Polar
plot_polar_uncertain(G_lon_q_nom, G_lon_q_unc_ARRAY, N, omega, ...
    'Polar: G_{lon,q}(j\omega) — \delta_{lon} \rightarrow q', ...
    'Nominal', 'Uncertain samples')

% --------------------------- Gtheta(s) ----------------------------------
% G_theta — Bode
plot_bode_unc(G_lon_theta_nom, G_lon_theta_unc_ARRAY, N, omega, ...
    'Bode: G_{lon,\theta}(j\omega) — \delta_{lon} \rightarrow \theta', ...
    '\omega [rad/s]', 'Magnitude [dB]', 'Nominal', 'Uncertain samples')

% G_theta — Polar
plot_polar_uncertain(G_lon_theta_nom, G_lon_theta_unc_ARRAY, N, omega, ...
    'Polar: G_{lon,\theta}(j\omega) — \delta_{lon} \rightarrow \theta', ...
    'Nominal', 'Uncertain samples')

%% Plotting functions

function plot_bode_unc(G_nom, G_unc_array, N, omega, plot_title, x_label, y_label, legend_nom, legend_unc)

figure

% --- Magnitude ---
ax1 = subplot(2,1,1); hold on
for i = 1:N
    [mag, ~] = bode(G_unc_array(:,:,i), omega);
    plot(omega, 20*log10(squeeze(mag)), ...   % ← plot normale
        'Color', [0.27 0.48 0.78 0.15], 'LineWidth', 0.8);
end
[mag_nom, ~] = bode(G_nom, omega);
h_nom = plot(omega, 20*log10(squeeze(mag_nom)), ...
    'Color', [0.85 0.15 0.15], 'LineWidth', 2.5);
set(ax1, 'XScale', 'log')   % ← imposta logaritmico esplicitamente
yline(0, 'k--', 'LineWidth', 1.2, 'Label', '0 dB', 'LabelHorizontalAlignment', 'left')
grid on
title(plot_title, 'FontSize', 13, 'FontWeight', 'bold')
xlabel(x_label,   'FontSize', 11)
ylabel(y_label,   'FontSize', 11)
h_patch = patch(NaN, NaN, [0.27 0.48 0.78], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend([h_patch, h_nom], ...
    {sprintf('%s (N=%d)', legend_unc, N), legend_nom}, ...
    'Location', 'southwest', 'FontSize', 10)

% --- Phase ---
ax2 = subplot(2,1,2); hold on
for i = 1:N
    [~, phase] = bode(G_unc_array(:,:,i), omega);
    plot(omega, squeeze(phase), ...   % ← plot normale
        'Color', [0.90 0.22 0.27 0.15], 'LineWidth', 0.8);
end
[~, phase_nom] = bode(G_nom, omega);
h_nom_ph = plot(omega, squeeze(phase_nom), ...
    'Color', [0.85 0.15 0.15], 'LineWidth', 2.5);
set(ax2, 'XScale', 'log')   % ← imposta logaritmico esplicitamente
yline(-180, 'k--', 'LineWidth', 1.2, 'Label', '-180°', 'LabelHorizontalAlignment', 'left')
yline( -90, 'k:',  'LineWidth', 1.0, 'Label',  '-90°', 'LabelHorizontalAlignment', 'left')
grid on
xlabel(x_label,      'FontSize', 11)
ylabel('Phase [deg]', 'FontSize', 11)
h_patch_ph = patch(NaN, NaN, [0.90 0.22 0.27], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend([h_patch_ph, h_nom_ph], ...
    {sprintf('%s (N=%d)', legend_unc, N), legend_nom}, ...
    'Location', 'southwest', 'FontSize', 10)
end


% Uncertain Polar Plot (with phase colormap)
function plot_polar_uncertain(G_nom, G_unc_array, N, omega, plot_title, legend_nom, legend_unc)

figure; hold on

for i = 1:N
    [mag, phase] = bode(G_unc_array(:,:,i), omega);
    re = squeeze(mag) .* cosd(squeeze(phase));
    im = squeeze(mag) .* sind(squeeze(phase));
    plot(re, im, 'Color', [0.27 0.48 0.78 0.15], 'LineWidth', 0.8)
end

[mag_nom, phase_nom] = bode(G_nom, omega);
re_nom = squeeze(mag_nom) .* cosd(squeeze(phase_nom));
im_nom = squeeze(mag_nom) .* sind(squeeze(phase_nom));
h_nom = plot(re_nom, im_nom, 'Color', [0.85 0.15 0.15], 'LineWidth', 2.5);

plot(-1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
xline(0, 'k-', 'LineWidth', 0.5, 'Alpha', 0.3)
yline(0, 'k-', 'LineWidth', 0.5, 'Alpha', 0.3)

grid on; axis equal
title(plot_title, 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Real',      'FontSize', 11)
ylabel('Imaginary', 'FontSize', 11)

h_patch = patch(NaN, NaN, [0.27 0.48 0.78], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend([h_patch, h_nom], ...
    {sprintf('%s (N=%d)', legend_unc, N), legend_nom}, ...
    'Location', 'best', 'FontSize', 10)
end

function uncertaintyRegionsPlot(G_nom, G_unc_array, N, omega_full, omega_regions, plot_title)
% UNCERTAINTYREGIONSPLOT  Polar plot with uncertainty disks at given frequencies
%
% Inputs:
%   G_nom         : nominal system (ss or tf)
%   G_unc_array   : uncertain samples array (from usample), size [:,:,N]
%   N             : number of samples
%   omega_full    : frequency vector for nominal curve [rad/s]
%   omega_regions : frequency points where to draw uncertainty disks [rad/s]
%   plot_title    : figure title string

N_regions = length(omega_regions);

% Evaluate frequency response at disk frequencies
G_nom_vals = squeeze(freqresp(G_nom, omega_regions));        % [N_regions x 1]
G_unc_vals = zeros(N_regions, N);
for i = 1:N
    G_unc_vals(:,i) = squeeze(freqresp(G_unc_array(:,:,i), omega_regions));
end

% Plot
figure; hold on; axis equal; grid on

% Nominal curve (full omega)
G_nom_full = squeeze(freqresp(G_nom, omega_full));
h_curve = plot(real(G_nom_full), imag(G_nom_full), ...
    'r-', 'LineWidth', 1.5);

% Uncertainty disks
for k = 1:N_regions
    G_nom_k = G_nom_vals(k);
    G_unc_k = G_unc_vals(k,:);

    % Cloud of uncertain points
    plot(real(G_unc_k), imag(G_unc_k), ...
        '.', 'Color', [0.27 0.48 0.78 0.4], 'MarkerSize', 8)

    % Disk radius = max distance from nominal
    r = max(abs(G_unc_k - G_nom_k));

    % Draw disk
    % start from nominal polar plot point and add points at 1-radius 
    % distance in 'all' direction to add points belonging to the circonference
    theta  = linspace(0, 2*pi, 100);
    h_circ = plot(real(G_nom_k) + r*cos(theta), ... 
        imag(G_nom_k) + r*sin(theta), ...
                  'b--', 'LineWidth', 1.2);

    % Nominal point
    h_nom_pt = plot(real(G_nom_k), imag(G_nom_k), ...
        'r*', 'MarkerSize', 8, 'LineWidth', 1.5);

    % Frequency label
    text(real(G_nom_k) + r*1.1, imag(G_nom_k), ...
        sprintf('\\omega=%.2f', omega_regions(k)), ...
        'FontSize', 8, 'Color', [0.2 0.2 0.2])
end

% Critical point
plot(-1, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2)
xline(0, 'k-', 'Alpha', 0.2)
yline(0, 'k-', 'Alpha', 0.2)

legend([h_curve, h_circ, h_nom_pt], ...
    {'Nominal G(j\omega)', 'Uncertainty disk', 'Nominal point'}, ...
    'Location', 'best', 'FontSize', 10)
xlabel('Real',      'FontSize', 11)
ylabel('Imaginary', 'FontSize', 11)
title(plot_title, 'FontSize', 13, 'FontWeight', 'bold')
end