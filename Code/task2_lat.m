%clc; clear; close all;

%% Load data
Aircraft = ANT_X_data();
g = Aircraft.g;

Yv_nom = Aircraft.Stability.Yv;
Yp_nom = Aircraft.Stability.Yp;
Lp_nom = Aircraft.Stability.Lp;
Yd_nom = Aircraft.Control.Yd;
Ld_nom = Aircraft.Control.Ld;

% Uncertainties: +- 3 sigma as deterministic bounds
% (Uncertainties are given in terms of standard deviations)
unc_Yv = 3 * Aircraft.Stability.Yv_sigma_pct;   
unc_Yp = 3 * Aircraft.Stability.Yp_sigma_pct;   
unc_Lp = 3 * Aircraft.Stability.Lp_sigma_pct;   
unc_Yd = 3 * Aircraft.Control.Yd_sigma_pct;    
unc_Ld = 3 * Aircraft.Control.Ld_sigma_pct;     

% Uncertain parameters
Yv_real = ureal('Yv_real', Yv_nom, 'Perc', unc_Yv);
Yp_real = ureal('Yp_real', Yp_nom, 'Perc', unc_Yp);
Lp_real = ureal('Lp_real', Lp_nom, 'Perc', unc_Lp);
Yd_real = ureal('Yd_real', Yd_nom, 'Perc', unc_Yd);
Ld_real = ureal('Ld_real', Ld_nom, 'Perc', unc_Ld);

fprintf('Uncertainty bounds \n')
fprintf('Yv: ±%.2f %%\n', unc_Yv)
fprintf('Yp: ±%.2f %%\n', unc_Yp)
fprintf('Lp: ±%.2f %%\n', unc_Lp)
fprintf('Yd: ±%.2f %%\n', unc_Yd)
fprintf('Ld: ±%.2f %%\n', unc_Ld)

%% Uncertain State-Space Model
% Uncertain Lateral Dynamics State Space Model
A_lat_unc = [Yv_real Yp_real g; 0 Lp_real 0; 0 1 0];
B_lat_unc = [Yd_real Ld_real 0]';
C_lat_unc = [0 1 0; 0 0 1];
D_lat_unc = [0 0]';

% Uncertain system
sys_lat_unc = ss(A_lat_unc, B_lat_unc, C_lat_unc, D_lat_unc);
sys_lat_unc.u = {'delta_lat'};
sys_lat_unc.y = {'p', 'phi'};

%% Nominal system
sys_lat_nom = sys_lat_unc.NominalValue;
sys_lat_nom.u = {'delta_lat'};
sys_lat_nom.y = {'p', 'phi'};

disp('Uncertain system'); 
disp(sys_lat_unc)

%%  Lateral Poles

% Nominal poles
disp('Lateral nominal poles')
[wn_lat, zeta_lat, p_lat] = damp(sys_lat_nom);
T = table(p_lat, wn_lat, zeta_lat, 'VariableNames', {'Pole','wn (rad/s)','zeta'});
disp(T)

% Uncertain poles (sampling) 
sys_lat_samples = usample(sys_lat_unc, 50);
poles_unc = pole(sys_lat_samples);
% disp('=== Uncertain poles (50 samples) ===')
% disp(poles_unc)

% delta_lat and p
figure('Name','Pole-Zero: p / delta_lat','NumberTitle','off');
pzplot(sys_lat_unc(1));
hold on
pzplot(sys_lat_nom(1));
grid on
legend('Uncertain','Nominal','Location','best')
title('Pole-Zero Map: \delta_{lat} \rightarrow p')

% delta_lat and phi
figure('Name','Pole-Zero: phi / delta_lat','NumberTitle','off');
pzplot(sys_lat_unc(2));
hold on
pzplot(sys_lat_nom(2));
grid on
legend('Uncertain','Nominal','Location','best')
title('Pole-Zero Map: \delta_{lat} \rightarrow \phi')

%% Lateral Zeros 

disp('Lateral nominal zeros: p/delta_lat')
disp(tzero(sys_lat_nom(1)))

disp('Lateral nominal zeros: phi/delta_lat')
disp(tzero(sys_lat_nom(2)))


%% Bode Plots

omega = logspace(-2, 3, 500);

% Bode: p / delta_lat
figure('Name', 'Bode: p / delta_lat', 'NumberTitle', 'off');
bode(sys_lat_unc(1), omega)
hold on
bode(sys_lat_nom(1), omega)
grid on
legend('Uncertain', 'Nominal', 'Location', 'northwest')
title('Bode: p / \delta_{lat}')

% Bode: phi / delta_lat
figure('Name', 'Bode: phi / delta_lat', 'NumberTitle', 'off');
bode(sys_lat_unc(2), omega)
hold on
bode(sys_lat_nom(2), omega)
grid on
legend('Uncertain', 'Nominal', 'Location', 'northwest')
title('Bode: \phi / \delta_{lat}')

%% Step response

figure('Name','Step Response','NumberTitle','off');
step(sys_lat_unc);
hold on
step(sys_lat_nom)
grid on
legend('Uncertain','Nominal','Location','northwest')
title('Lateral step Response – Uncertain and Nominal Models')

