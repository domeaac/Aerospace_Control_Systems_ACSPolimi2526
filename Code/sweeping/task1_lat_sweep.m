clc; clear; close all
PathsInit

Aircraft  = ANT_X_data();

Yv = Aircraft.Stability.Yv;
Yp = Aircraft.Stability.Yp;
Lp = Aircraft.Stability.Lp;
Yd = Aircraft.Control.Yd;
Ld = Aircraft.Control.Ld;

A_nom = [Yv Yp Aircraft.g; 0 Lp 0; 0 1 0];
B_nom = [Yd Ld 0]';
C_nom = [0 1 0; 0 0 1];
D_nom = [0 0]';

sys_lat_nom = ss(A_nom, B_nom, C_nom, D_nom, ...
    'StateName',  {'v (m/s)','p (rad/s)','phi (rad)'}, ...
    'InputName',  {'delta_lat'}, ...
    'OutputName', {'p','phi'});

omega_n_nom = 10;
xi_nom      = 0.9;

s = tf('s');
F_2_nom = omega_n_nom^2 / (s^2 + 2*xi_nom*omega_n_nom*s + omega_n_nom^2);
S_2_nom = 1 - F_2_nom;

[pk_nom, ~]       = getPeakGain(S_2_nom);
[~, ~, ~, wc_nom] = margin(S_2_nom);

fprintf('Nominal references:  wc = %.3f rad/s,  pk(S) = %.4f\n', wc_nom, pk_nom);

M_vec       = pk_nom .* (1:0.05:1.5);
omega_b_vec = wc_nom .* (0.8:0.05:1.3);
alpha_vec   = deg2rad(20) ./ (0.05:0.0025:0.1);

REQ_T_set      = 0.47;
REQ_delta      = 0.05;
REQ_gamma      = 1.0;
REQ_overshoot  = stepinfo(F_2_nom).Overshoot;

nM  = numel(M_vec);
nWb = numel(omega_b_vec);
nA  = numel(alpha_vec);

T_set_grid     = NaN(nA, nWb, nM);
maxDelta_grid  = NaN(nA, nWb, nM);
gamma_grid     = NaN(nA, nWb, nM);
overshoot_grid = NaN(nA, nWb, nM);
pk_S_grid      = NaN(nA, nWb, nM);

total = nA * nWb * nM;
fprintf('Starting sweep: %d x %d x %d = %d designs\n', nA, nWb, nM, total);
fprintf('Estimated time: %.0f - %.0f min\n\n', total*5/60, total*20/60);

t_start = tic;

for ia = 1:nA
    for iw = 1:nWb
        for im = 1:nM
            fprintf('  [%3d/%3d] alpha=%.3f, wb=%5.2f, M=%.2f  ...  ', ...
                (ia-1)*nWb*nM + (iw-1)*nM + im, total, ...
                alpha_vec(ia), omega_b_vec(iw), M_vec(im));

            [T_set_grid(ia,iw,im), maxDelta_grid(ia,iw,im), gamma_grid(ia,iw,im), ...
             overshoot_grid(ia,iw,im), pk_S_grid(ia,iw,im)] = ...
                task1_lat_func(M_vec(im), omega_b_vec(iw), alpha_vec(ia), sys_lat_nom);

            fprintf('Ts=%.3fs  dmax=%.4frad  gam=%.3f  OS=%.2f%%  pkS=%.3f\n', ...
                T_set_grid(ia,iw,im), maxDelta_grid(ia,iw,im), gamma_grid(ia,iw,im), ...
                overshoot_grid(ia,iw,im), pk_S_grid(ia,iw,im));
        end
    end
end

fprintf('\nSweep done in %.1f min\n', toc(t_start)/60);

save('task1_lat_sweep_results.mat', ...
    'M_vec','omega_b_vec','alpha_vec', ...
    'T_set_grid','maxDelta_grid','gamma_grid','overshoot_grid','pk_S_grid');

T_set_raw    = T_set_grid;
maxDelta_raw = maxDelta_grid;
gamma_raw    = gamma_grid;

gamma_mask = gamma_grid > REQ_gamma | ~isfinite(gamma_grid);
T_set_grid(gamma_mask)     = NaN;
maxDelta_grid(gamma_mask)  = NaN;
gamma_grid(gamma_mask)     = NaN;
overshoot_grid(gamma_mask) = NaN;
pk_S_grid(gamma_mask)      = NaN;

[alpha_asc, ia_sort] = sort(alpha_vec);
T_set_grid     = T_set_grid(ia_sort, :, :);
maxDelta_grid  = maxDelta_grid(ia_sort, :, :);
gamma_grid     = gamma_grid(ia_sort, :, :);
overshoot_grid = overshoot_grid(ia_sort, :, :);
pk_S_grid      = pk_S_grid(ia_sort, :, :);

metrics    = {T_set_grid,    maxDelta_grid,                 gamma_grid,  overshoot_grid,  pk_S_grid                  };
fig_names  = {'Settling Time','Max Actuator Deflection',    'H-inf Norm','Overshoot',     'S Peak'                   };
cb_labs    = {'T_{set} [s]', '|\delta_{lat}|_{max} [rad]', '\gamma',    'Overshoot [%]', '||S||_\infty [-]'         };
reqs       = {REQ_T_set,     REQ_delta,                    REQ_gamma,   REQ_overshoot,   Inf                        };
cmaps      = {'parula',      'hot',                        'cool',      'cool',          'cool'                     };
txt_colors = {[0 1 0],      [0 1 1],                      [1 1 0],     [0 1 1],         [1 0.6 0]                  };

M_nom     = pk_nom;
alpha_nom = deg2rad(20) / 0.05;
pick_ticks = @(v) v(unique(round(linspace(1, numel(v), min(numel(v), 10)))));

fsolve_opts = optimoptions('fsolve', 'Display', 'off');
xi_of_M = @(M) fsolve(@(xi) M - norm(tf([1, 2*xi*omega_n_nom, 0], ...
    [1, 2*xi*omega_n_nom, omega_n_nom^2]), inf), 0.4, fsolve_opts);
xi_vec = arrayfun(xi_of_M, M_vec);
fmt_M_tick = @(m) sprintf('%.3g (\\xi=%.2f)', m, xi_vec(M_vec == m));

[IA, IW, IM] = ndgrid(1:nA, 1:nWb, 1:nM);
xf_M  = M_vec(IM(:));
xf_wb = omega_b_vec(IW(:));
xf_a  = alpha_asc(IA(:));

for m = 1:numel(metrics)
    zf  = metrics{m}(:);
    ok  = isfinite(zf);
    bad = ok & (zf > reqs{m});

    figure('Name', fig_names{m}, 'NumberTitle','off', ...
           'Position', [80+(m-1)*40, 80+(m-1)*40, 960, 700]);

    if ~any(ok)
        title([fig_names{m} ' - no valid designs']);
        continue
    end

    clims = [min(zf(ok)) max(zf(ok))];
    if clims(1) == clims(2); clims = clims + [-0.1 0.1]; end

    hold on;
    scatter3(xf_M(ok), xf_wb(ok), xf_a(ok), 150, zf(ok), 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);
    colormap(cmaps{m}); clim(clims);
    cb = colorbar(); cb.Label.String = cb_labs{m};

    if any(bad)
        scatter3(xf_M(bad), xf_wb(bad), xf_a(bad), 250, 'x', ...
            'LineWidth',2.5, 'MarkerEdgeColor','k');
    end

    for k = find(ok)'
        text(xf_M(k), xf_wb(k), xf_a(k), sprintf('%.3g', zf(k)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',7, 'FontWeight','bold', 'Color',txt_colors{m});
    end

    if M_nom      >= min(M_vec)       && M_nom      <= max(M_vec)       && ...
       wc_nom     >= min(omega_b_vec) && wc_nom     <= max(omega_b_vec) && ...
       alpha_nom  >= min(alpha_asc)   && alpha_nom  <= max(alpha_asc)
        scatter3(M_nom, wc_nom, alpha_nom, 600, 'p', 'filled', ...
            'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth',1.5);
    end
    hold off;

    xtv = pick_ticks(M_vec); ytv = pick_ticks(omega_b_vec); ztv = pick_ticks(alpha_asc);
    set(gca, 'XTick',xtv, 'YTick',ytv, 'ZTick',ztv);
    xticklabels(arrayfun(fmt_M_tick, xtv, 'UniformOutput',false));
    yticklabels(arrayfun(@(v)sprintf('%.3g',v), ytv, 'UniformOutput',false));
    zticklabels(arrayfun(@(v)sprintf('%.3g',v), ztv, 'UniformOutput',false));

    xlabel('M_{lat}  (\xi)',  'FontSize',11);
    ylabel('\omega_b [rad/s]','FontSize',11);
    zlabel('\alpha_{lat}',    'FontSize',11);
    grid on; box on; view(45, 25);
    title(sprintf('%s\n color=value  |  X=fails req(%.3g)  |  *=nominal', ...
          fig_names{m}, reqs{m}), 'FontSize',11);
end

T_set_raw         = T_set_raw(ia_sort, :, :);
maxDelta_raw      = maxDelta_raw(ia_sort, :, :);
gamma_raw         = gamma_raw(ia_sort, :, :);
overshoot_raw     = overshoot_grid;

any_valid  = isfinite(T_set_raw) | isfinite(maxDelta_raw) | isfinite(gamma_raw);
f_gamma    = gamma_raw    > REQ_gamma     | ~isfinite(gamma_raw);
f_Tset     = T_set_raw    > REQ_T_set     | ~isfinite(T_set_raw);
f_delta    = maxDelta_raw > REQ_delta     | ~isfinite(maxDelta_raw);
f_overshoot = overshoot_raw > REQ_overshoot | ~isfinite(overshoot_raw);
n_fail     = double(f_gamma) + double(f_Tset) + double(f_delta) + double(f_overshoot);

n_fail_flat = n_fail(:);
valid_flat  = any_valid(:);

figure('Name','Requirements Satisfaction','NumberTitle','off', ...
       'Position',[200, 200, 960, 700]);
hold on;

colors_pass = [0.18 0.72 0.23];
colors_fail = [0.85 0.10 0.10];

pass = valid_flat & (n_fail_flat == 0);
fail = valid_flat & (n_fail_flat >  0);

if any(fail)
    scatter3(xf_M(fail), xf_wb(fail), xf_a(fail), 200, ...
        'filled', 'MarkerFaceColor',colors_fail, 'MarkerEdgeColor','k', 'LineWidth',0.5);
    for k = find(fail)'
        text(xf_M(k), xf_wb(k), xf_a(k), sprintf('%d fail', n_fail_flat(k)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',7, 'FontWeight','bold', 'Color',[1 0.5 0]);
    end
end
if any(pass)
    scatter3(xf_M(pass), xf_wb(pass), xf_a(pass), 200, ...
        'filled', 'MarkerFaceColor',colors_pass, 'MarkerEdgeColor','k', 'LineWidth',0.5);
    for k = find(pass)'
        text(xf_M(k), xf_wb(k), xf_a(k), 'OK', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',7, 'FontWeight','bold', 'Color',[0 0.4 0]);
    end
end

if M_nom     >= min(M_vec)       && M_nom     <= max(M_vec)       && ...
   wc_nom    >= min(omega_b_vec) && wc_nom    <= max(omega_b_vec) && ...
   alpha_nom >= min(alpha_asc)   && alpha_nom <= max(alpha_asc)
    scatter3(M_nom, wc_nom, alpha_nom, 600, 'p', 'filled', ...
        'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth',1.5);
end
hold off;

xtv = pick_ticks(M_vec); ytv = pick_ticks(omega_b_vec); ztv = pick_ticks(alpha_asc);
set(gca, 'XTick',xtv, 'YTick',ytv, 'ZTick',ztv);
xticklabels(arrayfun(fmt_M_tick, xtv, 'UniformOutput',false));
yticklabels(arrayfun(@(v)sprintf('%.3g',v), ytv, 'UniformOutput',false));
zticklabels(arrayfun(@(v)sprintf('%.3g',v), ztv, 'UniformOutput',false));

xlabel('M_{lat}  (\xi)',  'FontSize',11);
ylabel('\omega_b [rad/s]','FontSize',11);
zlabel('\alpha_{lat}',    'FontSize',11);
grid on; box on; view(45, 25);
title(sprintf('All Requirements  (T_{set}\\leq%.2fs, |\\delta|\\leq%.3frad, \\gamma\\leq%.1f)\n green=ALL OK  |  red>=1 fail  |  label=# failures  |  *=nominal', ...
      REQ_T_set, REQ_delta, REQ_gamma), 'FontSize',11);
