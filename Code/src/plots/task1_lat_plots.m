% plot chunks

if plot_nd
    fprintf(' \n Transfer function G_p(s) : p / \delta_{lat}')
    disp(G_p_nom)
    zpk(G_p_nom)

    fprintf(' \n Transfer function G_phi(s) : phi / \delta_{lat}')
    disp(G_phi_nom)
    zpk(G_phi_nom)

    figure % Plot Zeros and Poles of nominal G_p(s) 
    pzplot(sys_lat_nom(1));
    grid on ; 
    title('G_p(s) pole-zero map');
    legend('Nominal model','location','best')

    figure % Plot Zeros and Poles of nominal G_phi(s) 
    pzplot(sys_lat_nom(2));
    grid on ; 
    title('G_{\phi}(s) pole-zero map');
    legend('Nominal model','location','best')
end


if plot_nd
    fprintf('G_p zeros: \n')
    disp(z_p)
    fprintf('G_p poles: \n\n')
    disp(p_p)
end

if plot_nd
    fprintf('G_phi zeros: \n')
    disp(z_phi)
    fprintf('G_phi poles: \n\n')
    disp(p_phi)
end

if plot_nd
    disp('Open-loop transmission zeros:');
    Table_zeros = table(zeros_lat, ...
        'VariableNames', {'Transmission Zeros'});
    disp(Table_zeros)

    disp('Open-loop poles:');
    Table_poles = table(poles_lat, ...
        'VariableNames', {'Poles'});
    disp(Table_poles)
    
    if all(real(poles_lat) < 0)
        disp('System is STABLE in open loop')
    else
        disp('System is UNSTABLE or marginally stable in open loop.')
    end
end

if plot_nd 
    figure('Name', 'Pole-Zero Map', 'NumberTitle', 'off');
    pzmap([G_p_nom; G_phi_nom]); % ricostruisco il vettore ottenuto facendo tf in task1_lat.m
    title('G Pole-Zero Map - Open Loop');
    grid on;
end

if plot_nd
    figure; 
    bode(sys_lat_nom)                                      
    grid on ; 
    legend('Nominal model','location','best')
end

if plot_nd
    figure
    step(sys_lat_nom)
    grid on ; 
    legend('Nominal model','location','best')
end

if plot_nd
    figure
    step(F_2_phi); % step response of the desired feedbacked system F
    grid on
end

if plot_nd
    figure
    bodemag(F_2_phi)
    hold on
    % Desired approx S(jw)
    bodemag(S_2_phi)
    grid on
    legend('F_2','S_2')
end

if plot_nd
    fprintf('A   (DC gain (w=0)):  %.4f = %.2f dB \n', dc_lat, 20*log10(abs(dc_lat)));
    fprintf('wc   (Crossover) : %.4f rad/s \n', wc_lat);
    fprintf('PG_lat   (Peak gain):       %.4f = %.2f dB at omega = %.4f rad/s\n', pk_lat, 20*log10(pk_lat), wpk_lat);
end

if plot_nd
    figure
    bodemag(1/W_q_lat); title('1/$W_{Q_{lat}}(j\omega)$', 'interpreter', 'latex');
    grid on
end

if plot_nd
    figure
    bodemag(F_phi)
    hold on
    bodemag(F_2_phi)
    title('F_{\phi}(j\omega)');
    legend('F_{\phi}', 'F_{2,\phi}', 'location', 'best')
    grid on
end

if plot_nd
    figure
    step(F_phi)
    hold on
    step(F_2_phi)
    grid on
    legend('hinfstruct','second order', 'location','southeast')
end

if plot_nd
    stepinfo(F_phi);
end

if plot_nd
    figure
    bodemag(S_lat)
    grid on
    hold on
    bodemag(1/W_p_lat)
    bodemag(S_2_phi)
    legend('S','1/W_p', 'S_2', 'location','southeast')
end

if plot_nd
    figure
    bodemag(Q_lat)
    grid on
    hold on
    bodemag(1/W_q_lat)
    legend('Q','1/W_q', 'location','southeast')
end

%if plot_nd
%    figure
%    margin(L_lat)
%    grid on
%end

if plot_nd
    figure
    lsim(Q_lat,u,t); % simulate \delta_{lat} = Q*u in time interval t
    hold on
    grid on
end