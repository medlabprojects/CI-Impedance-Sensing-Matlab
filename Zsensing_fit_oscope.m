%%

% import data collected with oscilloscope
load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\testing\F1e1_3+4_04-03-2018.mat')

t_pulse = Z_pulse(:,1);

% smooth in two stages -> before/after jump
smooth_span = 0.08;
ind_jump = find(t_pulse>2.36e-6,1); % time picked by looking at Z_all
Z_pulse_smooth = Z_pulse(:,2);
Z_pulse_smooth(1:ind_jump)   = smooth(t_pulse(1:ind_jump),   Z_pulse(1:ind_jump,2),  8, 'rloess');
Z_pulse_smooth(ind_jump:end) = smooth(t_pulse(ind_jump:end), Z_pulse(ind_jump:end,2), smooth_span, 'rloess');

% resample at a fixed time interval (also ensures monotonically increases)
t_step = 0.5e-6; % [s]
t_interp_start = 8e-6; % [s]
t_interp = t_interp_start:t_step:max(t_pulse);
Z_pulse_interp = interp1(t_pulse, Z_pulse_smooth, t_interp);


%%
%-------------------------------------------------%
%  Compute estimates for Rp, Cp, & Ro (Method 2)  %
%-------------------------------------------------%

Zdot  = diff(Z_pulse_interp) ./ diff(t_interp);
ln_Zodot = log(Zdot);
Zodot_linfit = polyfit(t_interp(1:end-1), ln_Zodot, 1);

% Cp = e^-b
cp = exp(-Zodot_linfit(2));

% Rp = -1/(Cp*m)
rp = -1/(Zodot_linfit(1) * cp);

% Ro = Zo - Zp
Zp = rp * (1 - exp(-t_interp./(rp*cp)));
ro_full = Z_pulse_interp - Zp;
ro_linfit = polyfit(t_interp, ro_full, 1);
ro = ro_linfit(2);

% compute Zo using estimated parameters
Zo_est = ro + rp*(1-exp(-t_pulse./(rp*cp)));
fit_error = Z_pulse_smooth - Zo_est;

% print out parameters
fprintf('\n Cp = %.2f nF \n Rp = %.1f ohms \n Ro = %.2f ohms \n', 1e9*cp, rp, ro);

%% Plot

hf = figure(6); clf(6);
subplotXmanyY_er(3,1);
hold on
plot(1e6*t_pulse,  Z_pulse(:,2), 'Color',[0 0 0 0.2], 'LineWidth', 0.5)
plot(1e6*t_pulse,  Z_pulse_smooth, ':c', 'LineWidth', 1.5)
plot(1e6*t_interp, Z_pulse_interp, 'm+')
plot(1e6*t_pulse, Zo_est, '--r')
% xlabel('Time [\mus]')
ylabel('Impedance [\Omega]', 'FontSize',12)
title('Z_m(t)')

subplotXmanyY_er(3,2);
hold on
plot(1e6*t_interp(1:end-1), ln_Zodot, '+m')
plot(1e6*t_interp, polyval(Zodot_linfit, t_interp), '--r')
% xlabel('Time [\mus]')
title('Linearized Fit')

subplotXmanyY_er(3,3);
plot(1e6*t_pulse, fit_error)
xlabel('Time [\mus]', 'FontSize',14)
ylabel('Impedance [\Omega]', 'FontSize',12)
title('\DeltaZ (Model - Z_m)')

set(hf,'Visible','on')