%clc; clear; close all;

%% Load data
Aircraft = ANT_X_data();
g = Aircraft.g;

Xu_nom = Aircraft.Stability.Xu;
Xq_nom = Aircraft.Stability.Xq;
Mq_nom = Aircraft.Stability.Mq;
Xd_nom = Aircraft.Control.Xd;
Md_nom = Aircraft.Control.Md;

% Uncertainties: +- 3 sigma as deterministic bounds
% (Uncertainties are given in terms of standard deviations)
unc_Xu = 3 * Aircraft.Stability.Xu_sigma_pct;   
unc_Xq = 3 * Aircraft.Stability.Xq_sigma_pct;   
unc_Mq = 3 * Aircraft.Stability.Mq_sigma_pct;   
unc_Xd = 3 * Aircraft.Control.Xd_sigma_pct;    
unc_Md = 3 * Aircraft.Control.Md_sigma_pct;     

% Uncertain parameters
Xu_real = ureal('Xu_real', Xu_nom, 'Perc', unc_Xu);
Xq_real = ureal('Xq_real', Xq_nom, 'Perc', unc_Xq);
Mq_real = ureal('Mq_real', Mq_nom, 'Perc', unc_Mq);
Xd_real = ureal('Xd_real', Xd_nom, 'Perc', unc_Xd);
Md_real = ureal('Md_real', Md_nom, 'Perc', unc_Md);


fprintf('Uncertainty bounds \n')
fprintf('Xu: ±%.2f %%\n', unc_Xu)
fprintf('Xq: ±%.2f %%\n', unc_Xq)
fprintf('Mq: ±%.2f %%\n', unc_Mq)
fprintf('Xd: ±%.2f %%\n', unc_Xd)
fprintf('Md: ±%.2f %%\n', unc_Md)

%% Uncertain State-Space Model
% Uncertain Longitudinal Dynamics State Space Model
A_lon_unc = [Xu_real Xq_real -g; 0 Mq_real 0; 0 1 0];
B_lon_unc = [Xd_real Md_real 0]';
C_lon_unc = [0 1 0; 0 0 1];
D_lon_unc = [0 0]';


% Uncertain system
sys_lon_unc = ss(A_lon_unc, B_lon_unc, C_lon_unc, D_lon_unc);
sys_lon_unc.u = {'delta_lon'};
sys_lon_unc.y = {'q', 'theta'};

%% Nominal system
sys_lon_nom = sys_lon_unc.NominalValue;
sys_lon_nom.u = {'delta_lon'};
sys_lon_nom.y = {'q', 'theta'};

disp('Uncertain system'); 
disp(sys_lon_unc)

%%  Longitudinal Poles

% Nominal poles
disp('Longitudinal nominal poles')
[wn_lon, zeta_lon, q_lon] = damp(sys_lon_nom);
T = table(q_lon, wn_lon, zeta_lon, 'VariableNames', {'Pole','wn (rad_s)','zeta'});
disp(T)

% Uncertain poles (sampling) 
sys_lon_samples = usample(sys_lon_unc, 50);
poles_unc = pole(sys_lon_samples);
% disp('=== Uncertain poles (50 samples) ===')
% disp(poles_unc)

% delta_lon and q
figure('Name','Pole-Zero: q / delta_lon','NumberTitle','off');
pzplot(sys_lon_unc(1));
hold on
pzplot(sys_lon_nom(1));
grid on
legend('Uncertain','Nominal','Location','best')
title('Pole-Zero Map: \delta_{lon} \rightarrow q')

% delta_lon and theta
figure('Name','Pole-Zero: theta / delta_lon','NumberTitle','off');
pzplot(sys_lon_unc(2));
hold on
pzplot(sys_lon_nom(2));
grid on
legend('Uncertain','Nominal','Location','best')
title('Pole-Zero Map: \delta_{lon} \rightarrow \theta')

%% Longitudinal Zeros 

disp('Longitudinal nominal zeros: q/delta_lon')
disp(tzero(sys_lon_nom(1)))

disp('Longitudinal nominal zeros: theta/delta_lon')
disp(tzero(sys_lon_nom(2)))


%% Bode Plots
omega = logspace(-2, 3, 500);

% Bode: q / delta_lon
figure('Name', 'Bode: q / delta_lon', 'NumberTitle', 'off');
bode(sys_lon_unc(1), omega)
hold on
bode(sys_lon_nom(1), omega)
grid on
legend('Uncertain', 'Nominal', 'Location', 'northwest')
title('Bode: q / \delta_{lon}')

% Bode: theta / delta_lon
figure('Name', 'Bode: theta / delta_lon', 'NumberTitle', 'off');
bode(sys_lon_unc(2), omega)
hold on
bode(sys_lon_nom(2), omega)
grid on
legend('Uncertain', 'Nominal', 'Location', 'northwest')
title('Bode: \theta / \delta_{lon}')

%% Step response — full time range
figure('Name','Step Response — Full','NumberTitle','off');
step(sys_lon_unc); hold on; step(sys_lon_nom);
grid on
legend('Uncertain','Nominal','Location','northwest')
title('Longitudinal Step Response – Full Time Range')

%% Step response — zoom 0-100 s  (log y para theta)
t_z     = linspace(0.01, 100, 500)';      % começa em 0.01 para evitar log(0)
y_nom_z = step(sys_lon_nom, t_z);         % [500 x 2] : col1=q  col2=theta

figure('Name','Step Response — Zoom 0-100 s','NumberTitle','off');
ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');

for k = 1:50
    yk = step(sys_lon_samples(:,:,k), t_z);
    plot   (ax1, t_z, yk(:,1),              'Color',[0.25 0.55 0.90 0.20]);
    semilogy(ax2, t_z, max(yk(:,2), 1e-10), 'Color',[0.25 0.55 0.90 0.20]);
end

hq_n = plot   (ax1, t_z, y_nom_z(:,1),              'r','LineWidth',1.5);
ht_n = semilogy(ax2, t_z, max(y_nom_z(:,2), 1e-10), 'r','LineWidth',1.5);

ylabel(ax1,'q [rad/s]');    xlabel(ax1,'Time (s)');
ylabel(ax2,'\theta [rad]'); xlabel(ax2,'Time (s)');
title(ax1,'q / \delta_{lon}   (0-100 s)');
title(ax2,'\theta / \delta_{lon}   (0-100 s  |  log scale)');

hq_u = plot   (ax1, NaN, NaN, 'Color',[0.25 0.55 0.90]);
ht_u = semilogy(ax2, NaN, NaN, 'Color',[0.25 0.55 0.90]);
legend(ax1, [hq_u hq_n], 'Uncertain','Nominal','Location','best');
legend(ax2, [ht_u ht_n], 'Uncertain','Nominal','Location','best');