%%
waveform = importLecroyWaveform('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2019-11-07\C4_z_mux_ad62000001.csv');
i1 = find(waveform.t > 7e-6, 1);
i2 = find(waveform.t > 1.49e-3, 1);

pulse.current = 100e-6; % [A] current


%% smooth and trim
t_offset = 5e-6;

smooth_span = 150; % [pts]
v_smooth = smooth(waveform.t, waveform.v, smooth_span);

pulse.t = waveform.t(i1:i2);
pulse.t = pulse.t - t_offset;
pulse.v = v_smooth(i1:i2);
pulse.z = pulse.v ./ pulse.current; % compute impedance

% resample at a fixed time interval (also ensures monotonically increases)
t_step = 10e-6; % [s]
t_interp_start = 10e-6; % [s]
% pulse.t_interp = t_interp_start:t_step:max(pulse.t);
pulse.t_interp = t_interp_start:t_step:450e-6;
pulse.z_interp = interp1(pulse.t, pulse.z, pulse.t_interp);




%%
%-------------------------------------------------%
%  Compute estimates for Rp, Cp, & Ro (Method 2)  %
%-------------------------------------------------%

Zdot  = diff(pulse.z_interp) ./ diff(pulse.t_interp);
ln_Zodot = log(Zdot);
Zodot_linfit = polyfit(pulse.t_interp(1:end-1), ln_Zodot, 1);

% Cp = e^-b
cp = exp(-Zodot_linfit(2));

% Rp = -1/(Cp*m)
rp = -1/(Zodot_linfit(1) * cp);

% Ro = Zo - Zp
Zp = rp * (1 - exp(-pulse.t_interp ./ (rp*cp)));
ro_full = pulse.z_interp - Zp;
ro_linfit = polyfit(pulse.t_interp, ro_full, 1);
ro = ro_linfit(2);

% compute Zo using estimated parameters
Zo_est = ro + rp*(1-exp(-pulse.t./(rp*cp)));
fit_error = pulse.z - Zo_est;

% print out parameters
fprintf('\n Cp = %.2f nF \n Rp = %.1f ohms \n Ro = %.2f ohms \n', 1e9*cp, rp, ro);



%% plot

figure(6); clf(6); hold on
plot(waveform.t*1e3, waveform.v, 'Color',[0 0 0 0.3])
plot(pulse.t*1e3, pulse.v, 'g', 'LineWidth', 1.5)

legend('Measured', 'Smoothed')
xlabel('Time [ms]')
ylabel('Voltage [V]')

%%
figure(7); clf(7); hold on
plot(pulse.t*1e3, pulse.z)
plot(pulse.t_interp*1e3, pulse.z_interp, '*')
plot(pulse.t*1e3, Zo_est, 'r')


xlabel('Time [ms]')
% ylabel('Voltage [V]')

%%
figure(8); clf(8); hold on
plot(pulse.t_interp(1:end-1)*1e3, ln_Zodot)
plot(pulse.t_interp(1:end-1)*1e3, polyval(Zodot_linfit, pulse.t_interp(1:end-1)))
