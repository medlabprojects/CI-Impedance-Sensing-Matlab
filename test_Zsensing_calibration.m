%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Verification of impedance sensing method (both calibration and measurement)
%
%   Falstad circuit simulation => http://tinyurl.com/y5cejrae
%
%
%   Trevor Bruns
%   October 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  'True' parameters

% ro_true = 5800;  % [ohms] open channel series resistance
% rp_true = 14500;  % [ohms] parallel resistance
% cp_true = 5.25e-9; % ~5.45nF

ro_true = 8750;     % [ohms] open channel series resistance
rp_true = 21500;    % [ohms] parallel resistance
cp_true = 1.621e-6; % [F]

t_start = 0;      % [s]
% t_end   = 750e-6; % [s] needs to be long enough that we reach steady-state
% t_pts   = 2000;   % number of time points
t_end   = 400e-3; % [s] needs to be long enough that we reach steady-state
t_pts   = 100e3;   % number of time points



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Calibration Using Open-Channel Impedance   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute 'true' response using model

t = linspace(t_start, t_end, t_pts)'; % [s] time points to evaluate at
Zo_true = ro_true + rp_true*(1-exp(-t./(rp_true*cp_true))); % [ohms] impedance in open channel


%%% Create simulated measurements with noise

t_step_sim = 8.5e-6; % [s] interval between measurements
% t_step_sim = 8.5e-6; % [s] interval between measurements
t_start_sim = 3e-6;  % [s] our first measurement is slightly after pulse begins
t_sim = [t_start_sim:t_step_sim:t_end]'; % [s] time points to evaluate at

sigma = 5; % [ohms] standard deviation of impedance measurements
noise = normrnd(0,sigma, size(t_sim));

Zo_sim = interp1(t, Zo_true, t_sim) + noise;


%-------------------------------------------------%
%  Compute estimates for Rp, Cp, & Ro (Method 1)  %
%-------------------------------------------------%

% linearized form: ln(Rt-Zo) = ln(Rp) + [-1/(Rp*Cp)]*t

% We first need to find the steady state resistance, Rt
rt = Zo_true(end); % for 'true' model we can just use the last point

n_avg = floor(20 / 100 * length(t_sim)); % average the last N% of points
ind_avg_sim = length(t_sim) - n_avg;
rt_sim = mean(Zo_sim(ind_avg_sim:end));


% truncate measurements to avoid taking log(negative #)
n_trunc = 90 / 100; % find index where Zo is N% of steady-state value
ind_trunc = find(Zo_true > (n_trunc *rt), 1);
t_trunc = t(1:ind_trunc);

ind_trunc_sim = find(Zo_sim > (n_trunc*rt_sim),1);
t_trunc_sim = t_sim(1:ind_trunc_sim);


% Next, compute the best-fit line of  the linearized data:
% ln(Rt-Zo) = ln(Rp) + -1/(Rp*Cp) * t 
%     Y     =   b    +     m      * X
Zo_lin = log(rt - Zo_true(1:ind_trunc));
lin_fit = [ones(ind_trunc,1), t_trunc] \ Zo_lin;
Zo_lin_fit = lin_fit(1) + lin_fit(2)*t_trunc;

Zo_lin_sim = log(rt_sim - Zo_sim(1:ind_trunc_sim));
lin_fit_sim = [ones(ind_trunc_sim,1), t_trunc_sim] \ Zo_lin_sim;
Zo_lin_fit_sim = lin_fit_sim(1) + lin_fit_sim(2)*t_trunc_sim;


% Rp = e^b
rp_cal = exp(lin_fit(1));

rp_cal_sim = exp(lin_fit_sim(1));


% Cp = -1/(Rp*m)
cp_cal = -1/(rp_cal * lin_fit(2));

cp_cal_sim = -1/(rp_cal_sim * lin_fit_sim(2));


% Ro = Rt - Rp
ro_cal = rt - rp_cal;
ro_cal_sim = rt_sim - rp_cal_sim;



%-------------------------------------------------%
%  Compute estimates for Rp, Cp, & Ro (Method 2)  %
%-------------------------------------------------%
Zodot_true    = diff(Zo_true(1:ind_trunc)) ./ diff(t_trunc);
ln_Zodot_true = log(Zodot_true);
Zodot_true_linfit = [ones(ind_trunc-1, 1), t_trunc(1:end-1)] \ ln_Zodot_true;


Zodot_sim  = diff(Zo_sim(1:ind_trunc_sim)) ./ diff(t_trunc_sim);
ln_Zodot_sim = log(Zodot_sim);
Zodot_sim_linfit = [ones(ind_trunc_sim-1, 1), t_trunc_sim(1:end-1)] \ ln_Zodot_sim;

% Cp = e^-b
cp_cal2     = exp(-Zodot_true_linfit(1));
cp_cal2_sim = exp(-Zodot_sim_linfit(1));

% Rp = -1/(Cp*m)
rp_cal2     = -1/(Zodot_true_linfit(2) * cp_cal2);
rp_cal2_sim = -1/(Zodot_sim_linfit(2)  * cp_cal2_sim);

% Ro = Zo - Zp
Zp_cal2     = rp_cal2 * (1 - exp(-t./(rp_cal2*cp_cal2)));
Zp_cal2_sim = rp_cal2_sim * (1 - exp(-t_sim./(rp_cal2_sim*cp_cal2_sim)));
ro_cal2_full     = Zo_true - Zp_cal2;
ro_cal2_sim_full = Zo_sim  - Zp_cal2_sim;
ro_cal2 = ro_cal2_full(1);
ro_cal2_sim = ro_cal2_sim_full(1);


% compute Zo using estimated parameters
Zo_cal      = ro_cal + rp_cal*(1-exp(-t./(rp_cal*cp_cal))); 
Zo_cal_sim  = ro_cal_sim + rp_cal_sim*(1-exp(-t./(rp_cal_sim*cp_cal_sim))); 
Zo_cal2     = ro_cal2 + rp_cal2*(1-exp(-t./(rp_cal2*cp_cal2))); 
Zo_cal2_sim = ro_cal2_sim + rp_cal2_sim*(1-exp(-t./(rp_cal2_sim*cp_cal2_sim))); 

% compute errors in the noisy fit params
ro_error  = abs(ro_true-ro_cal_sim);
rp_error  = abs(rp_true-rp_cal_sim);
cp_error  = abs(cp_true-cp_cal_sim);
ro_error2 = abs(ro_true-ro_cal2_sim);
rp_error2 = abs(rp_true-rp_cal2_sim);
cp_error2 = abs(cp_true-cp_cal2_sim);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Closed-Channel Impedance Measurement (Zm)   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulse_width = 100e-6; % [s] length of positive portion of the constant-current pulse used for measurement
dt = -0.75e-6;    % [s] error between measured (assumed) time and true time
rc_true = 1500; % [ohms] true closed-channel resistance

% trim down time arrays (remove first 2.5us to avoid 'negative' times due to dt)
t_pulse_range = find(t>=2.5e-6, 1) : find(t>=pulse_width, 1);
t_pulse_sim_range = find(t_sim>=2.5e-6, 1) : find(t_sim>=pulse_width, 1);
t_pulse = t(t_pulse_range);
t_pulse_sim = t_sim(t_pulse_sim_range);

Zm_true = rc_true + Zo_true(t_pulse_range);             % 'true' measured closed-channel impedance
Zm_true_dt = rc_true + interp1(t, Zo_true, t_pulse-dt); % 'true' measured closed-channel impedance with timing uncertainty
Zm_sim = rc_true + Zo_sim(t_pulse_sim_range);
Zm_sim_dt  = rc_true + interp1(t, Zo_true, t_pulse_sim-dt) + noise(t_pulse_sim_range); % simulated closed-channel impedance measurements with noise and timing uncertainty (dt)

% compute the difference between the closed and open channel impedances
dZ_true    = Zm_true    - Zo_true(t_pulse_range);
dZ_true_dt = Zm_true_dt - Zo_true(t_pulse_range);
dZ_sim     = Zm_sim - interp1(t, Zo_cal_sim,  t_pulse_sim);
dZ_sim_dt  = Zm_sim_dt - interp1(t, Zo_cal, t_pulse_sim);

% figure; plot(t_pulse_sim, dZ_sim_dt); hold on; plot(t_pulse_sim, smooth(t_pulse_sim, dZ_sim_dt, 5, 'rlowess'),'--r')
% dZ_sim_dt = smooth(t_pulse_sim, dZ_sim_dt, 5, 'rlowess');


% if assuming no timing error, Rc = mean(dZ)
rc_sim = mean(dZ_sim);
rc_true_dt_uncorrected = mean(dZ_true_dt);
rc_sim_dt_uncorrected  = mean(dZ_sim_dt);


%%% determine the timing uncertainty (dt)
% first, compute the derivative of dZ wrt time
dZdot_true    = diff(dZ_true) ./ diff(t_pulse);
dZdot_true_dt = diff(dZ_true_dt) ./ diff(t_pulse);
dZdot_sim_dt  = diff(dZ_sim_dt) ./ diff(t_pulse_sim);

% if negative then dt is negative, but must flip regardless to ensure no instances of log(negative #)
if sum(dZdot_true_dt) < 0
    dt_sign_true_dt = -1;
else
    dt_sign_true_dt = 1;
end
dZdot_true_dt = abs(dZdot_true_dt);

if sum(dZdot_sim_dt) < 0
    dt_sign_sim_dt = -1;
else
    dt_sign_sim_dt = 1;
end
dZdot_sim_dt = abs(dZdot_sim_dt);


% next, linearize and find least-squares solution
% ln(dZdot*Cp) = ln[e^(dt/(Rp*Cp)) - 1]  +  -1/(Rp*Cp) * t
%      Y       =        b                +      m      * X
Y_Zm_linfit     = log(dZdot_true_dt*cp_cal);
Y_Zm_linfit_sim = log(dZdot_sim_dt*cp_cal_sim);

linfit_Zm_true_dt = polyfit( t_pulse(1:end-1),     Y_Zm_linfit,      1);
linfit_Zm_sim_dt  = polyfit( t_pulse_sim(1:end-1), Y_Zm_linfit_sim , 1);

% dt = Rp*Cp*ln(e^b+1)
dt_est_Zm_true_dt = dt_sign_true_dt * rp_cal*cp_cal*log(exp(linfit_Zm_true_dt(2))+1);
dt_est_Zm_sim_dt  = dt_sign_sim_dt  * rp_cal_sim*cp_cal_sim*log(exp(linfit_Zm_sim_dt(2))+1);

% Rc = Rp*e^(-t/(Rp*Cp)) - Rp*e^((dt-t)/(Rp*Cp)) - dZ
% should have zero slope, so just take the mean
rc_true_dt_corrected   = mean( dZ_true_dt + rp_cal*exp((dt_est_Zm_true_dt-t_pulse)./(rp_cal*cp_cal)) - rp_cal*exp(-t_pulse./(rp_cal*cp_cal)) );
rc_sim_dt_corrected    = mean( dZ_sim_dt + rp_cal_sim*exp((dt_est_Zm_sim_dt-t_pulse_sim)./(rp_cal_sim*cp_cal_sim)) - rp_cal_sim*exp(-t_pulse_sim./(rp_cal_sim*cp_cal_sim)) );


% create fitted lines using the parameter estimates found
linfit_dZ_true_dt = rc_true_dt_corrected - rp_cal*exp((dt_est_Zm_true_dt - t_pulse)./(rp_cal*cp_cal)) + rp_cal*exp(-t_pulse./(rp_cal*cp_cal));
linfit_dZ_sim_dt = rc_sim_dt_corrected - rp_cal_sim*exp((dt_est_Zm_sim_dt - t_pulse_sim)./(rp_cal_sim*cp_cal_sim)) + rp_cal_sim*exp(-t_pulse_sim./(rp_cal_sim*cp_cal_sim));


%% Create tables of parameters

Ro = [ro_true; ro_cal;  ro_cal_sim;  ro_error;  100*ro_error/ro_true;...
               ro_cal2; ro_cal2_sim; ro_error2; 100*ro_error2/ro_true];

Rp = [rp_true; rp_cal;  rp_cal_sim;  rp_error;  100*rp_error/rp_true;...
               rp_cal2; rp_cal2_sim; rp_error2; 100*rp_error2/rp_true];

Cp = [1e9*cp_true; 1e9*cp_cal;  1e9*cp_cal_sim;  1e9*cp_error;  100*cp_error/cp_true;...
                   1e9*cp_cal2; 1e9*cp_cal2_sim; 1e9*cp_error2; 100*cp_error2/cp_true]; % [nF]

method_names = {'True'; 'Cal1 (True)'; 'Cal1 (Simulated)'; 'Cal1 Sim Error (abs)'; 'Cal1 Sim Error (%)';...
                        'Cal2 (True)'; 'Cal2 (Simulated)'; 'Cal2 Sim Error (abs)'; 'Cal2 Sim Error (%)'};

param_table = table(Ro, Rp, Cp, 'RowNames',method_names);


% Rc
Rc = [rc_true; rc_sim; rc_true_dt_uncorrected; rc_sim_dt_uncorrected; rc_true_dt_corrected; rc_sim_dt_corrected];
Rc_err_abs = abs(Rc - rc_true);
Rc_err_pct = 100 * Rc_err_abs ./ rc_true;
Rc_table = table(Rc, Rc_err_abs, Rc_err_pct, 'RowNames', {'True'; 'Sim'; 'True (w/dt)'; 'Sim  (w/dt)'; 'True (dt-corrected)'; 'Sim  (dt-corrected)'});



%%
%%%%%%%%%%%%%%%
%%%  Plot   %%%
%%%%%%%%%%%%%%%

if exist('h_fig','var')
    if isvalid(h_fig)
        close(h_fig)
    end
end

h_fig = figure('name','Impedance Sensing Simulation','Visible','off');

% Impedance vs Time
h_ax(1) = subplot_er(3,2,[1,3]);
hold on
plot(1e6*t, 1e-3*Zo_true, 'k', 'LineWidth',1.5)
plot(1e6*t, 1e-3*Zo_cal, '--g', 'LineWidth',1.2)
plot(1e6*t_sim, 1e-3*Zo_sim, 'b*')
plot(1e6*t, 1e-3*Zo_cal_sim, '--r', 'LineWidth',1.3)
plot(1e6*t_sim(ind_avg_sim),   1e-3*Zo_sim(ind_avg_sim), 'ms', 'MarkerSize',20)
plot(1e6*t_sim(ind_trunc_sim), 1e-3*Zo_sim(ind_trunc_sim), 'ms', 'MarkerSize',20)
title(sprintf('Open-Channel Impedance vs Time  (t_{step,sim} = %.1f \\mus, \\sigma = %.1f \\Omega)', 1e6*t_step_sim, sigma), 'FontSize',13)
ylabel('Z_o [k\Omega]')
ylim([0 h_ax(1).YLim(2)])
legend('Z_{true}', 'Z_{cal,true}', 'Z_{sim}', 'Z_{cal,sim}', 'Location','SE')


% Linearized Impedance vs Time
h_ax(2) = subplot_er(3,2,5);
hold on
plot(1e6*t_trunc, Zo_lin, 'k', 'LineWidth',1.5)
plot(1e6*t_trunc, Zo_lin_fit, '--g', 'LineWidth',1.2)
plot(1e6*t_trunc_sim, Zo_lin_sim, 'b*')
plot(1e6*t_trunc_sim, Zo_lin_fit_sim, '--r', 'LineWidth',1.3)
legend('True', 'Calibrated (True)', 'Simulated', 'Calibrated (Simulated)', 'Location','SE')
title('Linearization to Estimate R_o, R_p, C_p')
xlabel('Time [\mus]')
ylabel('$ln(Z_{o,ss} - Z_o)$','Interpreter','latex', 'FontSize',12)


% Closed-channel measurement
h_ax(3) = subplot_er(3,2,2); 
hold on
plot(1e6*t,       1e-3*Zo_true,   'k', 'LineWidth',2)
plot(1e6*t,       1e-3*Zo_cal_sim, ':b', 'LineWidth',1.5)
plot(1e6*t_pulse, 1e-3*Zm_true,   ':g', 'LineWidth',1.5)
plot(1e6*t_pulse, 1e-3*Zm_true_dt,':r', 'LineWidth',1.5)
plot(1e6*t_pulse_sim, 1e-3*Zm_sim_dt, 'm*')
xlim([0,100])
ylim([0, h_ax(3).YLim(2)])
ylabel('k\Omega')
legend('Z_{o,true}', 'Z_{o,cal,sim}', 'Z_{m,true}', 'Z_{m,true,dt}', 'Z_{m,sim}', 'Location','SE')
title('Determining R_c from Closed-Channel Impedance Measurement (Z_m)', 'FontSize',13)

% dZ
h_ax(4) = subplot_er(3,2,4); 
hold on
plot(1e6*t_pulse,     dZ_true,   'k', 'LineWidth',2.5)
plot(1e6*t_pulse_sim, dZ_sim,     'go', 'MarkerFaceColor','g')
plot(1e6*t_pulse,     dZ_true_dt, '.b', 'LineWidth',2)
plot(1e6*t_pulse_sim, dZ_sim_dt,  'ms', 'MarkerFaceColor','m')
plot(1e6*t_pulse, linfit_dZ_true_dt, '--c', 'LineWidth',2) 
plot(1e6*t_pulse_sim, linfit_dZ_sim_dt, ':m', 'LineWidth',2)
plot(1e6*[0;t_pulse_sim], repmat(rc_sim, [1,length(t_pulse_sim)+1]), '--g', 'LineWidth',2)
plot(1e6*[0;t_pulse], repmat(rc_true_dt_uncorrected, [1,length(t_pulse)+1]), ':c', 'LineWidth',2)
plot(1e6*[0;t_pulse_sim], repmat(rc_sim_dt_uncorrected, [1,length(t_pulse_sim)+1]), '--m', 'LineWidth',2)
xlim([0,100])
% xlabel('Time [\mus]')
ylabel('\DeltaZ = Z_m - Z_o [\Omega]')
legend('\DeltaZ_{true}', '\DeltaZ_{sim}', '\DeltaZ_{true,dt}', '\DeltaZ_{sim,dt}', 'Location','NE')
title(sprintf('Correcting for Sampling Time Offset (\\Deltat = %.2f \\mus)', 1e6*dt))


% linearized 
h_ax(5) = subplot_er(3,2,6); 
hold on
plot(0, log(abs(exp(dt./(rp_true*cp_true))-1)), 'ko', 'MarkerFaceColor','k')
plot(1e6*t_pulse(1:end-1),     Y_Zm_linfit,     'b.')
plot(1e6*t_pulse_sim(1:end-1), Y_Zm_linfit_sim, 'ms', 'MarkerFaceColor','m')
plot(1e6*[0;t_pulse],     polyval(linfit_Zm_true_dt, [0;t_pulse]),      '--c', 'LineWidth',2)
plot(1e6*[0;t_pulse_sim], polyval(linfit_Zm_sim_dt,  [0;t_pulse_sim]),  '--m', 'LineWidth',2)

xlim([0,100])
xlabel('Time [\mus]', 'FontSize',13)
ylabel('$ln(C_{p}\dot{\Delta Z})$','Interpreter','latex', 'FontSize',12)
legend('Correct \Deltat Intercept', 'True (with \Deltat)', 'Simulated (with \Deltat)', 'Location','NE')
title('Linearization to Estimate \Deltat')



% Parameter Tables
% Ro, Rp, Cp
TString = evalc('format bank; disp(param_table); format shortG;'); % Get the table in string form.
TString = strrep(TString,'<strong>','\bf'); % Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName'); % Get a fixed-width font.

% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth, 'Units','Normalized', 'Position',[0.15 0.6 0.23 0.18]);


% Rc
TString = evalc('format bank; disp(Rc_table); format shortG;'); % Get the table in string form.
TString = strrep(TString,'<strong>','\bf'); % Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName'); % Get a fixed-width font.

% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth, 'Units','Normalized', 'Position',[0.68 0.71 0.25 0.12]);

h_fig.WindowState = 'maximized';
set(h_fig,'Visible','on')
linkaxes(h_ax(1:2), 'x');
linkaxes(h_ax(3:4), 'x');


%%

% figure(2);clf(2);
% subplot(2,1,1)
% hold on
% n = 30;
% plot(t_sim(1:n-1), diff(Zo_sim(1:n)), 'k', 'LineWidth',2)
% xlimits = get(gca,'XLim');
% line(xlimits, [0 0])
% % plot(t, ro_cal2_full, 'k', 'LineWidth',2)
% % plot(t_sim, ro_cal2_sim_full, 'r*');
% 
% subplot(2,1,2)
% hold on
% plot(t_trunc(1:end-1), ln_Zodot_true, 'k', 'LineWidth',2)
% plot(t_trunc_sim(1:end-1), ln_Zodot_sim, 'r*');
% xlim(xlimits)