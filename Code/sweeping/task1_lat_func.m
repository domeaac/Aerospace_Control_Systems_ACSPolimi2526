function [T_set, max_delta, gam, overshoot, pk_S] = task1_lat_func(M_lat, omega_b_lat, alpha_lat, sys_lat_nom)

T_set = NaN;  max_delta = NaN;  gam = NaN;  overshoot = NaN;  pk_S = NaN;

try
    s = tf('s');

    R_phi = tunablePID('R_phi', 'P');
    R_phi.u = 'e_phi';
    R_phi.y = 'p_0';

    R_p = tunablePID2('R_p', 'PID');
    R_p.c.Value  = 0;     R_p.c.Free  = false;
    R_p.b.Value  = 1;     R_p.b.Free  = false;
    R_p.Tf.Value = 0.02;  R_p.Tf.Free = false;
    R_p.u = {'p_0', 'p'};
    R_p.y = {'delta_lat'};

    SumOuter = sumblk('e_phi = phi_0 - phi');

    A_s   = 1e-3;
    w_tau = 50;

    W_p_inv = makeweight(A_s, omega_b_lat, M_lat, 0, 1);
    W_p = 1 / W_p_inv;
    W_p.u = 'e_phi';
    W_p.y = 'z_1';

    W_q = alpha_lat * (s + w_tau*1e-3) / (s + w_tau);
    W_q.u = 'delta_lat';
    W_q.y = 'z_2';

    CL0 = connect(R_p, R_phi, sys_lat_nom, SumOuter, W_p, W_q, ...
        'phi_0', {'z_1','z_2'}, {'delta_lat','e_phi','phi','p_0','p'});

    opt = hinfstructOptions('Display', 'off', 'RandomStart', 5);
    [CL, gam, ~] = hinfstruct(CL0, opt);

    F_phi = getIOTransfer(CL, 'phi_0', 'phi');
    Q_lat = getIOTransfer(CL, 'phi_0', 'delta_lat');

    info      = stepinfo(F_phi);
    T_set     = info.SettlingTime;
    overshoot = info.Overshoot;

    t = linspace(0, 6, 1e4);
    u = deg2rad( 10*(t > 1 & t <= 3) - 10*(t > 3 & t <= 5) );
    max_delta = max(abs(lsim(Q_lat, u, t)));

    S_cl  = getIOTransfer(CL, 'phi_0', 'e_phi');
    pk_S  = getPeakGain(S_cl);

catch ME
    warning('task1_lat_func: failed [M=%.2f wb=%.2f a=%.3f]: %s', ...
        M_lat, omega_b_lat, alpha_lat, ME.message);
end
end
