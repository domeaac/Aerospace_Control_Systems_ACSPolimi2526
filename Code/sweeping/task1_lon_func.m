function [T_set, max_delta, gam, overshoot, pk_S] = task1_lon_func(M_lon, omega_b_lon, alpha_lon, sys_lon_nom)

T_set = NaN;  max_delta = NaN;  gam = NaN;  overshoot = NaN;  pk_S = NaN;

try
    s = tf('s');

    R_theta = tunablePID('R_theta', 'P');
    R_theta.u = 'e_theta';
    R_theta.y = 'q_0';

    R_q = tunablePID2('R_q', 'PID');
    R_q.c.Value  = 0;     R_q.c.Free  = false;
    R_q.b.Value  = 1;     R_q.b.Free  = false;
    R_q.Tf.Value = 0.02;  R_q.Tf.Free = false;
    R_q.u = {'q_0', 'q'};
    R_q.y = {'delta_lon'};

    SumOuter = sumblk('e_theta = theta_0 - theta');

    A_s   = 1e-3;
    w_tau = 50;

    W_p_inv = makeweight(A_s, omega_b_lon, M_lon, 0, 1);
    W_p = 1 / W_p_inv;
    W_p.u = 'e_theta';
    W_p.y = 'z_1';

    W_q = alpha_lon * (s + w_tau*1e-3) / (s + w_tau);
    W_q.u = 'delta_lon';
    W_q.y = 'z_2';

    CL0 = connect(R_q, R_theta, sys_lon_nom, SumOuter, W_p, W_q, ...
        'theta_0', {'z_1','z_2'}, {'delta_lon','e_theta','theta','q_0','q'});

    opt = hinfstructOptions('Display', 'off', 'RandomStart', 5);
    [CL, gam, ~] = hinfstruct(CL0, opt);

    F_theta = getIOTransfer(CL, 'theta_0', 'theta');
    Q_lon   = getIOTransfer(CL, 'theta_0', 'delta_lon');

    info      = stepinfo(F_theta, 'SettlingTimeThreshold', 0.02);
    T_set     = info.SettlingTime;
    overshoot = info.Overshoot;

    t = linspace(0, 6, 1e4);
    u = deg2rad( 10*(t > 1 & t <= 3) - 10*(t > 3 & t <= 5) );
    max_delta = max(abs(lsim(Q_lon, u, t)));

    S_cl  = getIOTransfer(CL, 'theta_0', 'e_theta');
    pk_S  = getPeakGain(S_cl);

catch ME
    warning('task1_lon_func: failed [M=%.2f wb=%.2f a=%.3f]: %s', ...
        M_lon, omega_b_lon, alpha_lon, ME.message);
end
end
