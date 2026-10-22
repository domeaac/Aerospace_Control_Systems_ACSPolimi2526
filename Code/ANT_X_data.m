function Aircraft = ANT_X_data()

    Aircraft.g = 9.81;  % m/s^2

    % Stability derivatives + uncertainty whose std (1-sigma) given as a percentage of the nominal value 
    Aircraft.Stability = struct( ...
        'Yv',           -0.1068, ...  % 1/s
        'Yv_sigma_pct',   4,  ...  % 
        'Yp',            -0.1192, ...  % m/(s·rad)
        'Yp_sigma_pct',   2,  ...  % 
        'Lp',            2.6478, ...  % 1/s
        'Lp_sigma_pct',   2,  ...  
        'Xu',           -0.1068, ...  % 1/s
        'Xu_sigma_pct',   4,     ...  % %
        'Xq',            0.1192, ...  % m/(s·rad)
        'Xq_sigma_pct',   2,     ...  % %
        'Mq',           -2.6478, ...  % 1/s
        'Mq_sigma_pct',   2      ...  % %
    );

    % Control derivatives + uncertainty whose std (1-sigma) given as a percentage of the nominal value 
    Aircraft.Control = struct( ...
        'Yd',           10.1647, ...  % m/s^2
        'Yd_sigma_pct',    1.5,  ...  % %
        'Ld',           450.7085, ...  % rad/s^2
        'Ld_sigma_pct',    1,  ...  % %
        'Xd',           -10.1647,  ...  % m/s^2
        'Xd_sigma_pct',    1.5,    ...  % %
        'Md',           450.7085,  ...  % rad/s^2
        'Md_sigma_pct',    1       ...  % %
    );

    % Actuator limits
    Aircraft.Limits = struct( ...
        'delta_lat_max',  1, ...  % normalised
        'delta_lat_min', -1  ...  % normalised
    );

end