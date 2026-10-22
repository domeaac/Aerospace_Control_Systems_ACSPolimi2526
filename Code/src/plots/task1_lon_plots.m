% plot chunks

if plot_nd
    fprintf(' \n Transfer function G_q(s) : q / \delta_{lon}')
    disp(G_q_nom)
    zpk(G_q_nom)

    fprintf(' \n Transfer function G_theta(s) : theta / \delta_{lon}')
    disp(G_theta_nom)
    zpk(G_theta_nom)

    figure % Plot Zeros and Poles of nominal G_q(s) 
    pzplot(sys_lon_nom(1));
    grid on ; 
    title('G_q(s) pole-zero map');
    legend('Nominal model','location','best')

    figure % Plot Zeros and Poles of nominal G_theta(s) 
    pzplot(sys_lon_nom(2));
    grid on ; 
    title('G_{\theta}(s) pole-zero map');
    legend('Nominal model','location','best')
end


if plot_nd
    fprintf('G_q zeros: \n')
    disp(z_q)
    fprintf('G_q poles: \n\n')
    disp(p_q)
end

if plot_nd
    fprintf('G_theta zeros: \n')
    disp(z_theta)
    fprintf('G_theta poles: \n\n')
    disp(p_theta)
end

if plot_nd
    disp('Open-loop transmission zeros:');
    Table_zeros = table(zeros_lon, ...
        'VariableNames', {'Transmission Zeros'});
    disp(Table_zeros)

    disp('Open-loop poles:');
    Table_poles = table(poles_lon, ...
        'VariableNames', {'Poles'});
    disp(Table_poles)
    
    if all(real(poles_lon) < 0)
        disp('System is STABLE in open loop')
    else
        disp('System is UNSTABLE or marginally stable in open loop.')
    end
end

if plot_nd
    figure('Name', 'Pole-Zero Map', 'NumberTitle', 'off');
    pzmap([G_q_nom; G_theta_nom]);
    title('G Pole-Zero Map - Open Loop');
    grid on;
end

if plot_nd
    figure; 
    bode(sys_lon_nom)                                      
    grid on ; 
    legend('Nominal model','location','best')
end

if plot_nd
    figure
    step(sys_lon_nom)
    grid on ;
    legend('Nominal model','location','best')
end

if plot_nd
    figure
    step(F_2_theta); % step response of the desired feedbacked system F
    grid on
end

if plot_nd
    figure
    bodemag(F_2_theta)
    hold on
    % Desired approx S(jw)
    bodemag(S_2_theta)
    grid on
    legend('F_2','S_2')
end

if plot_nd
    fprintf('A   (DC gain (w=0)):  %.4f = %.2f dB \n', dc_lon, 20*log10(abs(dc_lon)));
    fprintf('wc   (Crossover) : %.4f rad/s \n', wc_lon);
    fprintf('PG_lon   (Peak gain):       %.4f = %.2f dB at omega = %.4f rad/s\n', pk_lon, 20*log10(pk_lon), wpk_lon);
end

if plot_nd
    figure
    bodemag(1/W_q_lon); title('1/$W_{Q_{lon}}(j\omega)$', 'interpreter', 'latex');
    grid on
end

if plot_nd
    figure
    bodemag(F_theta)
    hold on
    bodemag(F_2_theta)
    title('F_{\theta}(j\omega)');
    legend('F_{\theta}', 'F_{2,\theta}', 'location', 'best')
    grid on
end

if plot_nd
    figure
    step(F_theta)
    hold on
    step(F_2_theta)
    grid on
    legend('hinfstruct','second order', 'location','southeast')
end

if plot_nd
    stepinfo(F_theta);
end

if plot_nd
    figure
    bodemag(S_lon)
    grid on
    hold on
    bodemag(1/W_p_lon)
    bodemag(S_2_theta)
    legend('S','1/W_p', 'S_2', 'location','southeast')
end

if plot_nd
    figure
    bodemag(Q_lon)
    grid on
    hold on
    bodemag(1/W_q_lon)
    legend('Q','1/W_q', 'location','southeast')
end

%if plot_nd
%    figure
%    margin(L_lon)
%    grid on
%end

if plot_nd
    figure
    lsim(Q_lon,u,t); % simulate \delta_{lon} = Q*u in time interval t
    hold on
    grid on
end