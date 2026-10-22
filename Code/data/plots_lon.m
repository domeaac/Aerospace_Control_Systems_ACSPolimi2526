%% task1_lon_sweep_plots.m
clc; close all;

load('task1_lon_sweep_results_30.mat');

%% ---- Parameters ----
omega_n_nom    = 10;
xi_nom         = 0.9;
s              = tf('s');
F_2_nom        = omega_n_nom^2 / (s^2 + 2*xi_nom*omega_n_nom*s + omega_n_nom^2);
S_2_nom        = 1 - F_2_nom;
[pk_nom, ~]    = getPeakGain(S_2_nom);
[~,~,~,wc_nom] = margin(S_2_nom);

REQ_T_set     = 0.47;
REQ_delta     = 0.05;
REQ_gamma     = 1.0;
REQ_overshoot = stepinfo(F_2_nom).Overshoot;

M_nom     = pk_nom;
alpha_nom = deg2rad(20) / 0.05;

%% ---- Sort alpha ascending ----
[alpha_asc, ia_sort] = sort(alpha_vec);
T_set_grid     = T_set_grid(ia_sort,:,:);
maxDelta_grid  = maxDelta_grid(ia_sort,:,:);
gamma_grid     = gamma_grid(ia_sort,:,:);
overshoot_grid = overshoot_grid(ia_sort,:,:);
pk_S_grid      = pk_S_grid(ia_sort,:,:);

%% ---- Raw copies before masking ----
T_set_raw    = T_set_grid;
maxDelta_raw = maxDelta_grid;
gamma_raw    = gamma_grid;

%% ---- Mask designs with gamma > 1 ----
gamma_mask = gamma_grid > REQ_gamma | ~isfinite(gamma_grid);
T_set_grid(gamma_mask)     = NaN;
maxDelta_grid(gamma_mask)  = NaN;
gamma_grid(gamma_mask)     = NaN;
overshoot_grid(gamma_mask) = NaN;
pk_S_grid(gamma_mask)      = NaN;

%% ---- Flat coordinate arrays ----
[IA, IW, IM] = ndgrid(1:numel(alpha_asc), 1:numel(omega_b_vec), 1:numel(M_vec));
xf_M  = M_vec(IM(:));
xf_wb = omega_b_vec(IW(:));
xf_a  = alpha_asc(IA(:));

%% ---- Helper functions ----
pick_ticks  = @(v) v(unique(round(linspace(1,numel(v),min(numel(v),10)))));
fsolve_opts = optimoptions('fsolve','Display','off');
xi_of_M     = @(M) fsolve(@(xi) M - norm(tf([1,2*xi*omega_n_nom,0], ...
                  [1,2*xi*omega_n_nom,omega_n_nom^2]),inf), 0.4, fsolve_opts);
xi_vec      = arrayfun(xi_of_M, M_vec);
fmt_M_tick  = @(m) sprintf('%.3g (\\xi=%.2f)', m, xi_vec(M_vec==m));

nominal_in_range = M_nom     >= min(M_vec)       && M_nom     <= max(M_vec)       && ...
                   wc_nom    >= min(omega_b_vec)  && wc_nom    <= max(omega_b_vec) && ...
                   alpha_nom >= min(alpha_asc)    && alpha_nom <= max(alpha_asc);

add_nominal = @(ax) scatter3(M_nom, wc_nom, alpha_nom, 600, 'p','filled', ...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5,'Parent',ax);

function set_ax_ticks(ax, M_vec, omega_b_vec, alpha_asc, fmt_M_tick)
    pick = @(v) v(unique(round(linspace(1,numel(v),min(numel(v),10)))));
    xtv = pick(M_vec); ytv = pick(omega_b_vec); ztv = pick(alpha_asc);
    set(ax,'XTick',xtv,'YTick',ytv,'ZTick',ztv);
    xticklabels(ax, arrayfun(fmt_M_tick, xtv,'UniformOutput',false));
    yticklabels(ax, arrayfun(@(v)sprintf('%.3g',v), ytv,'UniformOutput',false));
    zticklabels(ax, arrayfun(@(v)sprintf('%.3g',v), ztv,'UniformOutput',false));
end

%% ================================================================
%% FIGURES 1-5 — 3D scatter per metric
%% ================================================================
metrics   = {T_set_grid,   maxDelta_grid,                 gamma_grid,  overshoot_grid,  pk_S_grid         };
fig_names = {'Settling Time','Max Actuator Deflection',   'H-inf Norm','Overshoot',     'S Peak'          };
cb_labs   = {'T_{set} [s]','|\delta_{lon}|_{max} [rad]', '\gamma',    'Overshoot [%]', '||S||_\infty'    };
reqs      = {REQ_T_set,    REQ_delta,                     REQ_gamma,   REQ_overshoot,   Inf               };
cmaps     = {'parula',     'hot',                         'cool',      'cool',          'cool'            };
txt_col   = {[0 .8 0],    [0 .9 .9],                     [1 1 0],     [0 .9 .9],       [1 .6 0]          };

for m = 1:numel(metrics)
    zf  = metrics{m}(:);
    ok  = isfinite(zf);
    bad = ok & (zf > reqs{m});

    figure('Name',fig_names{m},'NumberTitle','off','Position',[60+m*25,60+m*25,960,700]);
    ax = gca; hold(ax,'on');

    if ~any(ok)
        title(ax,[fig_names{m} ' — no valid designs']); continue
    end

    clims = [min(zf(ok)) max(zf(ok))];
    if clims(1)==clims(2), clims = clims+[-0.1 0.1]; end

    scatter3(xf_M(ok),xf_wb(ok),xf_a(ok), 150, zf(ok),'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.5,'Parent',ax);
    colormap(ax,cmaps{m}); clim(ax,clims);
    cb = colorbar(ax); cb.Label.String = cb_labs{m};

    if any(bad)
        scatter3(xf_M(bad),xf_wb(bad),xf_a(bad), 250,'x', ...
            'LineWidth',2.5,'MarkerEdgeColor','k','Parent',ax);
    end

    for k = find(ok)'
        text(xf_M(k),xf_wb(k),xf_a(k), sprintf('%.3g',zf(k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7,'FontWeight','bold','Color',txt_col{m},'Parent',ax);
    end

    if nominal_in_range, add_nominal(ax); end
    hold(ax,'off');

    set_ax_ticks(ax, M_vec, omega_b_vec, alpha_asc, fmt_M_tick);
    xlabel(ax,'M_{lon} (\xi)','FontSize',14);
    ylabel(ax,'\omega_b [rad/s]','FontSize',14);
    zlabel(ax,'\alpha_{lon}','FontSize',14);
    grid(ax,'on'); box(ax,'on'); view(ax,45,25);
    title(ax,sprintf('%s\n color=value  |  X=violates req(%.3g)  |  %s=nominal', ...
        fig_names{m},reqs{m},char(9733)),'FontSize',14);
end

%% ================================================================
%% FIGURE 6 — Requirements satisfaction map
%% ================================================================
any_valid   = isfinite(T_set_raw)|isfinite(maxDelta_raw)|isfinite(gamma_raw);
f_gamma     = gamma_raw    > REQ_gamma     | ~isfinite(gamma_raw);
f_Tset      = T_set_raw    > REQ_T_set     | ~isfinite(T_set_raw);
f_delta     = maxDelta_raw > REQ_delta     | ~isfinite(maxDelta_raw);
f_overshoot = overshoot_grid>REQ_overshoot | ~isfinite(overshoot_grid);
n_fail      = double(f_gamma)+double(f_Tset)+double(f_delta)+double(f_overshoot);

n_fail_flat = n_fail(:);
valid_flat  = any_valid(:);
pass = valid_flat & (n_fail_flat==0);
fail = valid_flat & (n_fail_flat>0);

figure('Name','Requirements Satisfaction','NumberTitle','off','Position',[250,250,960,700]);
ax2 = gca; hold(ax2,'on');

if any(fail)
    scatter3(xf_M(fail),xf_wb(fail),xf_a(fail), 200,'filled', ...
        'MarkerFaceColor',[0.85 0.10 0.10],'MarkerEdgeColor','k','LineWidth',0.5,'Parent',ax2);
    for k = find(fail)'
        text(xf_M(k),xf_wb(k),xf_a(k),sprintf('%d fail',n_fail_flat(k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7,'FontWeight','bold','Color',[1 0.5 0],'Parent',ax2);
    end
end
if any(pass)
    scatter3(xf_M(pass),xf_wb(pass),xf_a(pass), 200,'filled', ...
        'MarkerFaceColor',[0.18 0.72 0.23],'MarkerEdgeColor','k','LineWidth',0.5,'Parent',ax2);
    for k = find(pass)'
        text(xf_M(k),xf_wb(k),xf_a(k),'OK', ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7,'FontWeight','bold','Color',[0 0.4 0],'Parent',ax2);
    end
end

if nominal_in_range, add_nominal(ax2); end
hold(ax2,'off');

set_ax_ticks(ax2, M_vec, omega_b_vec, alpha_asc, fmt_M_tick);
xlabel(ax2,'M_{lon} (\xi)','FontSize',14);
ylabel(ax2,'\omega_b [rad/s]','FontSize',14);
zlabel(ax2,'\alpha_{lon}','FontSize',14);
grid(ax2,'on'); box(ax2,'on'); view(ax2,45,25);
title(ax2,sprintf( ...
    'All Requirements  (T_{set}\\leq%.2fs  |  |\\delta|\\leq%.3frad  |  \\gamma\\leq%.1f)\n green=ALL OK  |  red=\\geq1 fail  |  label=#failures  |  %s=nominal', ...
    REQ_T_set,REQ_delta,REQ_gamma,char(9733)),'FontSize',14);

%% ================================================================
%% FIGURE 7 — 2D heatmap slices at most feasible alpha
%% ================================================================
pass_3d     = reshape(pass,[numel(alpha_asc),numel(omega_b_vec),numel(M_vec)]);
[~,best_ia] = max(sum(pass_3d,[2 3]));

figure('Name','2D slice — best alpha','NumberTitle','off','Position',[300,100,1200,800]);
sub_metrics = {T_set_grid, maxDelta_grid,              gamma_grid,  overshoot_grid};
sub_names   = {'T_{set} [s]','|\delta|_{max} [rad]',   '\gamma',    'Overshoot [%]'};
sub_reqs    = {REQ_T_set,    REQ_delta,                 REQ_gamma,   REQ_overshoot };
sub_cmaps   = {'parula',     'hot',                     'cool',      'cool'        };

for m = 1:4
    ax3 = subplot(2,2,m);
    slice = squeeze(sub_metrics{m}(best_ia,:,:));   % [nWb x nM]
    imagesc(ax3, M_vec, omega_b_vec, slice);
    set(ax3,'YDir','normal');
    colormap(ax3,sub_cmaps{m}); colorbar(ax3);
    hold(ax3,'on');
    if isfinite(sub_reqs{m})
        contour(ax3,M_vec,omega_b_vec,slice,[sub_reqs{m} sub_reqs{m}], ...
            'r-','LineWidth',2,'DisplayName',sprintf('req=%.3g',sub_reqs{m}));
    end
    plot(ax3,M_nom,wc_nom,'wp','MarkerSize',14,'MarkerFaceColor','w', ...
         'MarkerEdgeColor','k','LineWidth',1.5);
    hold(ax3,'off');
    xlabel(ax3,'M_{lon}','FontSize',13);
    ylabel(ax3,'\omega_b [rad/s]','FontSize',13);
    title(ax3,sprintf('%s  (\\alpha=%.3g)',sub_names{m},alpha_asc(best_ia)),'FontSize',13);
    grid(ax3,'on');
end
sgtitle(sprintf('2D slices at \\alpha_{lon} = %.4f (most feasible)',alpha_asc(best_ia)),'FontSize',14);

%% ================================================================
%% FPRINTF — Summary
%% ================================================================
fprintf('\n========================================\n')
fprintf('           SWEEP SUMMARY\n')
fprintf('========================================\n')
fprintf('  Total designs evaluated : %d\n', numel(T_set_raw))
fprintf('  gamma <= 1 (feasible)   : %d\n', sum(~gamma_mask(:)))
fprintf('  ALL requirements met    : %d\n', sum(pass))
fprintf('\n  Requirements:\n')
fprintf('    T_set     <= %.3f s\n',   REQ_T_set)
fprintf('    |delta|   <= %.4f rad\n', REQ_delta)
fprintf('    gamma     <= %.2f\n',     REQ_gamma)
fprintf('    Overshoot <= %.2f %%\n',  REQ_overshoot)

if any(pass)
    pass_idx     = find(pass);
    [~,best_rel] = min(gamma_raw(pass_idx));
    best_abs     = pass_idx(best_rel);
    fprintf('\n  Best design (min gamma among passing):\n')
    fprintf('    M_lon      = %.4f  (xi=%.4f)\n', xf_M(best_abs),  xi_of_M(xf_M(best_abs)))
    fprintf('    omega_b    = %.4f rad/s\n',       xf_wb(best_abs))
    fprintf('    alpha      = %.4f\n',             xf_a(best_abs))
    fprintf('    T_set      = %.4f s\n',   T_set_raw(best_abs))
    fprintf('    |delta|max = %.4f rad\n', maxDelta_raw(best_abs))
    fprintf('    gamma      = %.4f\n',     gamma_raw(best_abs))
    fprintf('    Overshoot  = %.4f %%\n',  overshoot_grid(best_abs))
else
    fprintf('\n  No design satisfies all requirements.\n')
    fprintf('  Consider relaxing REQ_T_set or REQ_delta.\n')
end
fprintf('========================================\n')