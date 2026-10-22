clc; clear; close all;

wn      = 10;
xi_vec  = linspace(0.1,1,15);
omega   = logspace(-1, 2, 800);
t       = linspace(0, 5, 2000);

s = tf('s');

cmap = lines(numel(xi_vec));

fig1 = figure('Name','Bode Magnitude: F and S','NumberTitle','off');
fig1.Position = [100 100 1100 500];

ax_F = subplot(1,2,1); hold on; grid on;
title('|F(j\omega)|  –  Complementary Sensitivity','FontSize',13)
xlabel('\omega  [rad/s]'); ylabel('Magnitude  [dB]')
xscale log; set(ax_F,'XLim',[omega(1) omega(end)])
yline(  0, 'k--', '0 dB', 'LabelHorizontalAlignment','left','FontSize',9)
yline(-3, 'k:' , '-3 dB','LabelHorizontalAlignment','left','FontSize',9)

ax_S = subplot(1,2,2); hold on; grid on;
title('|S(j\omega)|  –  Sensitivity','FontSize',13)
xlabel('\omega  [rad/s]'); ylabel('Magnitude  [dB]')
xscale log; set(ax_S,'XLim',[omega(1) omega(end)])
yline(  0, 'k--', '0 dB', 'LabelHorizontalAlignment','left','FontSize',9)

leg_str = cell(1, numel(xi_vec));

for k = 1:numel(xi_vec)
    xi = xi_vec(k);

    F = wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    S = 1 - F;

    [magF] = squeeze(bode(F, omega));
    [magS] = squeeze(bode(S, omega));

    plot(ax_F, omega, 20*log10(magF), 'Color', cmap(k,:), 'LineWidth', 1.8)
    plot(ax_S, omega, 20*log10(magS), 'Color', cmap(k,:), 'LineWidth', 1.8)

    leg_str{k} = sprintf('\\xi = %.3g', xi);
end

legend(ax_F, leg_str, 'Location','southwest','FontSize',10)
legend(ax_S, leg_str, 'Location','southeast','FontSize',10)

xline(ax_F, wn, 'k-.', '\omega_n', 'LabelVerticalAlignment','bottom','FontSize',9)
xline(ax_S, wn, 'k-.', '\omega_n', 'LabelVerticalAlignment','bottom','FontSize',9)

fig2 = figure('Name','Step Response: F and S','NumberTitle','off');
fig2.Position = [100 640 1100 420];

ax_Fstep = subplot(1,2,1); hold on; grid on;
title('Step Response of F(s)','FontSize',13)
xlabel('Time  [s]'); ylabel('Amplitude  [-]')
yline(1, 'k--', 'steady state','LabelHorizontalAlignment','left','FontSize',9)

ax_Sstep = subplot(1,2,2); hold on; grid on;
title('Step Response of S(s)','FontSize',13)
xlabel('Time  [s]'); ylabel('Amplitude  [-]')
yline(0, 'k--','','LabelHorizontalAlignment','left','FontSize',9)

for k = 1:numel(xi_vec)
    xi = xi_vec(k);

    F = wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    S = 1 - F;

    yF = step(F, t);
    yS = step(S, t);

    plot(ax_Fstep, t, yF, 'Color', cmap(k,:), 'LineWidth', 1.8)
    plot(ax_Sstep, t, yS, 'Color', cmap(k,:), 'LineWidth', 1.8)
end

legend(ax_Fstep, leg_str, 'Location','southeast','FontSize',10)
legend(ax_Sstep, leg_str, 'Location','northeast','FontSize',10)

xi_fine = linspace(0.05, 1.5, 300);
Ms_vec  = zeros(1, numel(xi_fine));

for k = 1:numel(xi_fine)
    xi = xi_fine(k);
    F  = wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    S  = 1 - F;
    Ms_vec(k) = getPeakGain(S);
end

figure('Name','Ms = ||S||_inf  vs  xi','NumberTitle','off');
fig3 = gcf; fig3.Position = [1220 100 540 380];
hold on; grid on;
plot(xi_fine, Ms_vec, 'b-', 'LineWidth', 2)
plot(xi_fine, 20*log10(Ms_vec), 'r--', 'LineWidth', 2)
yline(1,    'k:',  'Ms = 1',    'FontSize', 9)
yline(1.4,  'k--', 'Ms = 1.4',  'FontSize', 9)
xline(0.5,  'b:',  '\xi = 0.5', 'FontSize', 9)
xline(0.707,'b--', '\xi = 1/\surd2','FontSize',9)

for k = 1:numel(xi_vec)
    xi  = xi_vec(k);
    F_k = wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    S_k = 1 - F_k;
    Ms_k = getPeakGain(S_k);
    plot(xi, Ms_k, 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor', 'k')
end

xlabel('\xi  (damping ratio)');
ylabel('Ms  [-]  (linear)  /  [dB]  (dashed)')
title('Sensitivity Peak  M_s = \|S\|_\infty  vs  \xi','FontSize',13)
legend([leg_str, {'Ms (linear)', 'Ms (dB)'}], ...
       'Location','northeast','FontSize',9)

fprintf('\n%s\n', repmat('=',1,60))
fprintf('  wn = %.4g rad/s\n', wn)
fprintf('%s\n', repmat('=',1,60))
fprintf('  %-8s  %-10s  %-12s  %-12s  %-12s\n', ...
        'xi', 'Ms=||S||inf', 'wb_F [rad/s]', 'OS_F [%]', 'Ts_F [s]')
fprintf('%s\n', repmat('-',1,60))

for k = 1:numel(xi_vec)
    xi = xi_vec(k);
    F  = wn^2 / (s^2 + 2*xi*wn*s + wn^2);
    S  = 1 - F;

    Ms   = getPeakGain(S);
    wb   = bandwidth(F);
    info = stepinfo(F);

    fprintf('  %-8.3g  %-10.4f  %-12.4f  %-12.4f  %-12.4f\n', ...
            xi, Ms, wb, info.Overshoot, info.SettlingTime)
end
fprintf('%s\n\n', repmat('=',1,60))
