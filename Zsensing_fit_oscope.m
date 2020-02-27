%% Import data collected with oscilloscope

%% Trial 1
% load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\testing\F1e1_3+4_04-03-2018.mat')
% 
% t_pulse = Z_pulse(:,1);
% 
% % smooth in two stages -> before/after jump
% smooth_span = 0.08;
% ind_jump = find(t_pulse>2.36e-6,1); % time picked by looking at Z_all
% ind_jump = 25;
% Z_pulse_smooth = Z_pulse(:,2);
% % Z_pulse_smooth(1:ind_jump)   = smooth(t_pulse(1:ind_jump),   Z_pulse(1:ind_jump,2),  8, 'rloess');
% Z_pulse_smooth(ind_jump:end) = smooth(t_pulse(ind_jump:end), Z_pulse(ind_jump:end,2), smooth_span, 'rloess');
% 
% % Resample at a fixed time interval (also ensures monotonically increases)
% t_step = 0.5e-6; % [s]
% t_interp_start = 8e-6; % [s]
% t_interp = t_pulse:t_step:max(t_pulse);
% Z_pulse_interp = interp1(t_pulse, Z_pulse_smooth, t_interp);


%% Trial 2
load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2019-10-31\z_mux_100us.mat')
smooth_span = 0.05;
ind_jump = 105;
t_pulse = Z_pulse(:,1);
Z_pulse_smooth = Z_pulse(:,2);
Z_pulse_smooth(1:ind_jump) = smoothdata(t_pulse(1:ind_jump));
Z_pulse_smooth(ind_jump:end) = smooth(t_pulse(ind_jump:end), Z_pulse(ind_jump:end,2), smooth_span, 'sgolay',2);

%% Resample at a fixed time interval (also ensures monotonically increases)
t_step = 6.5e-6; % [s]
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
% Z vs t
ax(1) = subplotXmanyY_er(3,1);
hold on
plot(1e6*t_pulse,  Z_pulse(:,2), 'b', 'LineWidth', 1.5)
% plot(1e6*t_pulse,  Z_pulse(:,2), 'Color',[0 0 0 0.2], 'LineWidth', 2.5)
% plot(1e6*t_pulse,  Z_pulse_smooth, ':c', 'LineWidth', 2.5)
plot(1e6*t_pulse, Zo_est, 'g', 'LineWidth', 2.5)
plot(1e6*t_interp, Z_pulse_interp, 'ro', 'MarkerSize',7, 'MarkerFaceColor','r')


ylim([-50 8000])
ylabel('\Omega', 'FontSize',12)
title('Measured Impedance', 'FontSize',14)
legend('Z_m(t)', 'Z_{fit}(t)', 'Measurement Points', 'FontSize',12, 'Location','se')


% Linearized
ax(2) = subplotXmanyY_er(3,2);
hold on
plot(1e6*t_interp(1:end-1), ln_Zodot, 'ro', 'MarkerSize',7, 'MarkerFaceColor','r')
plot(1e6*t_interp, polyval(Zodot_linfit, t_interp), 'g', 'LineWidth', 2)


title('Linearized Fit', 'FontSize',14)


% Model Error
ax(3) = subplotXmanyY_er(3,3);
idx = find(t_pulse > t_interp(1), 1);
plot(1e6*t_pulse(idx:end), fit_error(idx:end), 'm', 'LineWidth',2)
xlabel('Time [\mus]', 'FontSize',12)
ylabel('\Delta\Omega', 'FontSize',12)
title('Residuals', 'FontSize',14)
legend('\DeltaZ (Z_m(t) - Z_{fit}(t)', 'FontSize',12, 'Location','se')

linkaxes(ax,'x')
xlim([-2 100])
set(hf,'Visible','on')



%% Plot for Dissertation

hf2 = figure(11); clf(11);

% current pulse
% ax2(1) = subplot_er(6,2,[1 3 5]);
ax2(1) = subplot_er(2,2,1);
hold on
I_pulse = [-20, -1e-9,  0, 99.99, 100, 104.99,  105, 199.99, 200, 220;
             0,   0   100,   100,   0,      0, -100,   -100,   0,   0];

plot(I_pulse(1,:),I_pulse(2,:), 'k', 'LineWidth',2)
yline(0,'k');

ylim([-120 120])
title('Current Pulse', 'FontSize',15)
ylabel('\muA', 'FontSize',14)
grid on


% voltage response
% ax2(2) = subplot_er(6,2,[7 9 11]);
ax2(2) = subplot_er(2,2,3);
hold on
idx = find(V_all(:,1) > V_pulse(end,1), 1);
plot(1e6*V_pulse(:,1),     1e3*V_pulse(:,2), 'Color',[0 0.4470 0.7410], 'LineWidth',2);
plot(1e6*V_all(idx:end,1), 1e3*V_all(idx:end,2), 'k', 'LineWidth',2);

yline(0,'k');

title('Measured Voltage Response', 'FontSize',14)
ylabel('mV', 'FontSize',12)
xlabel('Time (\mus)', 'FontSize',12)
ylim([-350 850])
linkaxes(ax2(1:2),'x')
xlim([-15 215])
grid on

% Z vs t
% ax2(3) = subplot_er(6,2,[2 4]);
ax2(3) = subplot_er(2,2,2);
hold on

hz(1) = plot(1e6*t_pulse,  Z_pulse(:,2),'Color',[0 0.4470 0.7410], 'LineWidth', 1.5);
% plot(1e6*t_pulse,  Z_pulse(:,2), 'Color',[0 0 0 0.2], 'LineWidth', 2.5)
% plot(1e6*t_pulse,  Z_pulse_smooth, ':c', 'LineWidth', 2.5)
hz(2) = plot(1e6*t_pulse, Zo_est, '--g', 'LineWidth', 1.5);
hz(3) = plot(1e6*t_interp, Z_pulse_interp, 'ro', 'MarkerSize',4, 'MarkerFaceColor','r');

grid on
xticks(0:25:100)
ylim([-50 8000])
ylabel('\Omega', 'FontSize',12)
title('Measured Impedance', 'FontSize',14)
legend(hz, {'Z_m(t)'; 'Z_{fit}(t)'; 'Measurements'}, 'FontSize',9.5, 'Location','se')


% Linearized
% ax2(4) = subplot_er(6,2,[6 8]);
ax2(4) = subplot_er(2,2,4);
hold on
plot(1e6*t_interp(1:end-1), ln_Zodot, 'ro', 'MarkerSize',4, 'MarkerFaceColor','r')
plot(1e6*t_interp, polyval(Zodot_linfit, t_interp), '--g', 'LineWidth', 1.5)
txt_params = {sprintf('R_a = %.0f \\Omega',ro); sprintf('R_p = %.0f \\Omega',rp); sprintf('C_p = %.2f nF',1e9*cp)};
annotation('textbox',[.77 .4 .02 .02],'String',txt_params,'FitBoxToText','on', 'FontSize',9.5, 'BackgroundColor','w');

ylabel('$\ln(\dot{Z})$', 'Interpreter','latex', 'FontSize',12)
grid on
title('Linearized Fit', 'FontSize',14)
xlabel('Time (\mus)', 'FontSize',12)



% % Model Error
% ax2(5) = subplot_er(6,2,[10 12]);
% idx = find(t_pulse > t_interp(1), 1);
% plot(1e6*t_pulse(idx:end), fit_error(idx:end), 'Color',[0.8500 0.3250 0.0980] , 'LineWidth',2)
% 
% grid on
% title('Residuals', 'FontSize',15)
% xlabel('Time [\mus]', 'FontSize',14)
% ylabel('\Delta\Omega', 'FontSize',14)
% legend('\DeltaZ = Z_m(t) - Z_{fit}(t)', 'FontSize',11, 'Location','se')

linkaxes(ax2(3:5),'x')
xlim([-2 100])
xticks([0:25:100])


set(hf2,'Visible','on')
