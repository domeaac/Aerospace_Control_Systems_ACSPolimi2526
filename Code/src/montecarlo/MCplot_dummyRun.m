%% Progressive Sigma Plot

sd_data = {sd_wc_phi,        sd_wc_theta, ...
           sd_xi_phi,        sd_xi_theta, ...
           sd_delta_lat_max, sd_delta_lon_max};

figure('Position', [100, 100, 1200, 800]);

titles = {'\omega_c^\phi [rad/s]', '\omega_c^\theta [rad/s]', ...
          '\xi_\phi [-]',          '\xi_\theta [-]', ...
          '\delta_{lat}^{max} [rad]', '\delta_{lon}^{max} [rad]'};

colors = {'#E63946','#457B9D','#2A9D8F','#E9C46A','#F4A261','#A8DADC'};

for j = 1:6
    subplot(2, 3, j);

    plot(1:k, sd_data{j}, 'Color', colors{j}, 'LineWidth', 1.5);

    % Converged mean as horizontal reference
    yline(sd_data{j}(end), 'k--', 'LineWidth', 1.2, ...
          'Label', sprintf('\\mu_{\\infty}=%.4f', sd_data{j}(end)), ...
          'LabelOrientation', 'horizontal');

    title(titles{j}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('N iterations'); ylabel('Progressive sigma \sigma');
    grid on; box on;
end

sgtitle(sprintf('Progressive stdev — CLT Convergence  (N = %d)', k), ...
        'FontSize', 14, 'FontWeight', 'bold');
% CLT: as i >, the sample mean converges to the mean of sigma (we should
% use k = infty).