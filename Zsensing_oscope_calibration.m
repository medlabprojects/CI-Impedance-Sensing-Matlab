%   Zsensing_oscope_calibration
%
%   Trevor Bruns
%   October 2019
%

%% user parameters

% filename = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Calibration testing\2019-10-23\F1_z_400ms00000.csv'; % E6-E7
% pulse.length  = 400e-3; % [s]
% pulse.current = 100e-6; % [A]

filename = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2019-10-29\C4_z_mux_500us00000_trimmed.csv';
pulse.length  = 500e-6; % [s]
pulse.current = 100e-6; % [A]

%% import and preprocess waveform captured by scope

% import
waveform = importLecroyWaveform(filename);

% smooth
smooth_span = 750; % [pts]
v_smooth = smooth(waveform.t, waveform.v, smooth_span);

% trim
trim_threshold = -1e-4; % [dV]
v_smooth_diff = diff(v_smooth);
i_pulseon = find(v_smooth_diff <= trim_threshold, 1); % find pulse start (voltage drops)
i_start = i_pulseon + find(v_smooth_diff(i_pulseon:end) > 0, 1); % find where voltage rise starts
i_pulseoff = length(v_smooth_diff) - find(flip(v_smooth_diff) >= -trim_threshold, 1); % find where pulse ends
i_end = find(waveform.t > (waveform.t(i_pulseoff) - 2e-6), 1); % trim back 2us to remove voltage fluctuations

pulse.t = waveform.t(i_start:i_end);
pulse.v = v_smooth(i_start:i_end);

% set pulse start as t = 0 seconds
pulse.t = pulse.t - waveform.t(i_pulseon);

% 'bias' up the voltage if <0
if pulse.v(1) < 0
    bias = pulse.v(1) - 0.1;
    pulse.v = pulse.v - bias;
else
    bias = 0;
end

% resample at a fixed time interval (also ensures monotonically increases)
t_step = 1.5e-6; % [s]
t_interp_start = 8e-6; % [s]
pulse.t_interp = t_interp_start:t_step:max(pulse.t);
pulse.z_interp = interp1(pulse.t, pulse.v, pulse.t_interp);

% compute impedance
pulse.z = pulse.v ./ pulse.current;



% % resample to reduce total number of data points
% npts = 5e3;
% resample.t = linspace(waveform.t(1), waveform.t(end), npts);
% resample.v = pchip(waveform.t, waveform.v, resample.t);
% 
% % smooth
% smooth_span = 10; % [pts]
% v_smooth = smooth(resample.t, resample.v, smooth_span);
% 
% % trim
% trim_threshold = 0.01; % [dV]
% v_smooth_diff = diff(v_smooth);
% i_start = find(v_smooth_diff >= trim_threshold, 1);
% i_end   = length(v_smooth_diff) - find(flip(v_smooth_diff) <= -trim_threshold, 1);
% pulse.t = resample.t(i_start:i_end)';
% pulse.v = v_smooth(i_start:i_end);
% 
% % set pulse start as t = 0 seconds
% pulse.t = pulse.t - pulse.t(1);
% 
% % compute impedance
% pulse.z = pulse.v ./ pulse.current;


%%
%-------------------------------------------------%
%  Compute estimates for Rp, Cp, & Ro (Method 1)  %
%-------------------------------------------------%

% linearized form: ln(Rt-Zo) = ln(Rp) + [-1/(Rp*Cp)]*t

% We first need to find the steady state resistance, Rt
n_avg = floor(20 / 100 * length(pulse.t)); % average the last N% of points
ind_avg = length(pulse.t) - n_avg;
rt = mean(pulse.z(ind_avg:end));


% truncate measurements to avoid taking log(negative #)
n_trunc = 90 / 100; % find index where Zo is N% of steady-state value
ind_trunc = find(pulse.z > (n_trunc*rt),1);
t_trunc = pulse.t(1:ind_trunc);


% Next, compute the best-fit line of  the linearized data:
% ln(Rt-Zo) = ln(Rp) + -1/(Rp*Cp) * t 
%     Y     =   b    +     m      * X
Zo_lin = log(rt - pulse.z(1:ind_trunc));
lin_fit = [ones(ind_trunc,1), t_trunc] \ Zo_lin;
Zo_lin_fit = lin_fit(1) + lin_fit(2)*t_trunc;


% Rp = e^b
rp_cal = exp(lin_fit(1))


% Cp = -1/(Rp*m)
cp_cal = -1/(rp_cal * lin_fit(2))


% Ro = Rt - Rp
ro_cal = rt - rp_cal

% compute Zo curve using estimated values
Zo_cal = ro_cal + rp_cal * (1 - exp(-pulse.t ./(rp_cal*cp_cal)));



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
% plot(resample.t*1e3, resample.v, 'r')
plot(waveform.t*1e3, v_smooth, 'g', 'LineWidth', 1.5)
% plot(resample.t*1e3, v_smooth, 'g', 'LineWidth', 1.5)

% legend('Measured', 'Calibrated')
xlabel('Time [ms]')
ylabel('Voltage [V]')

%%
figure(7); clf(7); hold on
plot(pulse.t*1e3, pulse.z/1e3, 'k', 'LineWidth',1.3)
% plot(pulse.t*1e3, Zo_cal/1e3, '--b', 'LineWidth', 1.2)

% legend('Measured', 'Calibrated')
xlabel('Time [ms]')
ylabel('Impedance [k\Omega]')