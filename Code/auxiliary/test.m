sys_lat_unc.u = {'delta_lat'};
sys_lat_unc.y = {'p','phi'};

CL_uncertain = connect(SumOuterLat, R_phi_tuned, R_p_tuned, sys_lat_unc, ...
                       {'phi_0'}, {'phi','p', 'delta_lat'});

opts = robOptions('Display', 'on', ...
                  'Sensitivity', 'on', ...
                  'VaryFrequency', 'on');

[stabmarg, wcu, info] = robstab(CL_uncertain, opts);
fprintf('System handles %.1f%% of modeled uncertainty.\n', stabmarg.LowerBound * 100);
