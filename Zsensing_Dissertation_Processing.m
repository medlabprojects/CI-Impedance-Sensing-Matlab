%% Import Impedance Data

directory_base = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data';

% trial 4 (for model characterization)
load(fullfile(directory_base, 'Zsense_2020-01-24_phantom_Flex24_5x20_EA-1-5_trial4.mat'));
trial4.A_raw = areas;
trial4.Z_raw = Z_sync;

% trial 7 (for model characterization)
load(fullfile(directory_base, 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial7.mat'));
trial7.A_raw = areas;
trial7.Z_raw = Z_sync;

% trial 6 (for testing model accuracy)
load(fullfile(directory_base, 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6.mat'));
trial6.A_raw = areas;
trial6.Z_raw = Z_sync;

% remove all other variables no longer needed
clear areas areas_raw idx_Z_sync pixels_per_meter polygons Z Z_sync Z_trial4 Z_trial6 Z_trial7

%%

directory_base = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data';
load(fullfile(directory_base, 'Zsense_dissertation-data_2-23.mat'));


%% Preprocessing
% remove bad impedance measurements (i.e. >100 Ohms below normal low-level baseline)

% manual and iterative process => end result is idx_inliers matrix for each trial
% ransac_max_distance = 0.8;
% ransac_iter = 2;
% ch = 4;
% [idx_inliers, coeffs, coeffs_lin] = ZsensingFitCurveRansac(trial7.Z(:,ch), trial7.A(:,ch), ransac_max_distance, ransac_iter, true);


% trial 6 used for validation, so use all values as inliers
trial6.idx_inliers = true(size(trial6.Z_raw));

% create cell arrays with inlier data
for ii=1:4
    trial4.A{ii,1} = trial4.A_raw( trial4.idx_inliers(:,ii), ii )';
    trial4.Z{ii,1} = trial4.Z_raw( trial4.idx_inliers(:,ii), ii )';

    trial7.A{ii,1} = trial7.A_raw( trial7.idx_inliers(:,ii), ii )';
    trial7.Z{ii,1} = trial7.Z_raw( trial7.idx_inliers(:,ii), ii )';

    trial6.A{ii,1} = trial6.A_raw( trial6.idx_inliers(:,ii), ii )';
    trial6.Z{ii,1} = trial6.Z_raw( trial6.idx_inliers(:,ii), ii )';
end


%% exponential moving average

alpha = 1/300;
EMA = @(x,last_ema) alpha*x + (1-alpha)*last_ema;

for i_ch=1:4
    trial4.Z_ema{i_ch,1} = zeros(1,numel(trial4.Z{i_ch}));
    trial4.Z_ema{i_ch}(1) = mean(trial4.Z{i_ch}(1:20)); % initialize using mean of first 20 samples

    trial7.Z_ema{i_ch,1} = zeros(1,numel(trial7.Z{i_ch}));
    trial7.Z_ema{i_ch}(1) = mean(trial7.Z{i_ch}(1:20)); % initialize using mean of first 20 samples

    trial6.Z_ema{i_ch,1} = zeros(1,numel(trial6.Z{i_ch}));
    trial6.Z_ema{i_ch}(1) = mean(trial6.Z{i_ch}(1:20)); % initialize using mean of first 20 samples

    for ii=2:numel(trial4.Z{i_ch})
        trial4.Z_ema{i_ch}(ii) = EMA(trial4.Z{i_ch}(ii), trial4.Z_ema{i_ch}(ii-1));
    end
    for ii=2:numel(trial7.Z{i_ch})
        trial7.Z_ema{i_ch}(ii) = EMA(trial7.Z{i_ch}(ii), trial7.Z_ema{i_ch}(ii-1));
    end
    for ii=2:numel(trial6.Z{i_ch})
        trial6.Z_ema{i_ch}(ii) = EMA(trial6.Z{i_ch}(ii), trial6.Z_ema{i_ch}(ii-1));
    end
end


%% current minimum

for i_ch=1:4
    for ii=1:numel(trial4.Z{i_ch})
        trial4.Z_min{i_ch,1}(ii) = min(trial4.Z{i_ch}(1:ii));
    end
    for ii=1:numel(trial7.Z{i_ch})
        trial7.Z_min{i_ch,1}(ii) = min(trial7.Z{i_ch}(1:ii));
    end
    for ii=1:numel(trial6.Z{i_ch})
        trial6.Z_min{i_ch,1}(ii) = min(trial6.Z{i_ch}(1:ii));
    end
end


%% determine bias value for each channel (i.e. open-channel impedance)

trial4.A_all = [];
trial4.Z_all = [];
trial7.A_all = [];
trial7.Z_all = [];
trial6.A_all = [];
trial6.Z_all = [];

% trial4.bias = [1950 1980 1995 1935];
% trial7.bias = [2300 2350 2350 2390];
for ii=1:4
    trial4.bias(ii) = min(trial4.Z{ii})-1;
    trial4.Z_bias{ii,1} = trial4.Z{ii} - trial4.bias(ii);
    trial4.Z_all = [trial4.Z_all, trial4.Z{ii}-trial4.bias(ii)];
    trial4.A_all = [trial4.A_all, trial4.A{ii}];

    trial7.bias(ii) = min(trial7.Z{ii})-1;
    trial7.Z_bias{ii,1} = trial7.Z{ii} - trial7.bias(ii);
    trial7.Z_all = [trial7.Z_all, trial7.Z{ii}-trial7.bias(ii)];
    trial7.A_all = [trial7.A_all, trial7.A{ii}];

    trial6.bias(ii) = min(trial6.Z_raw(:,ii))-1;
    trial6.Z_bias{ii,1} = trial6.Z{ii} - trial6.bias(ii);
    trial6.Z_all = [trial6.Z_all, trial6.Z_raw(:,ii)'-trial6.bias(ii)];
    trial6.A_all = [trial6.A_all, trial6.A_raw(:,ii)'];
end

t4Z1 = trial4.Z{1} - trial4.bias(1);
t4Z2 = trial4.Z{2} - trial4.bias(2);
t4Z3 = trial4.Z{3} - trial4.bias(3);
t4Z4 = trial4.Z{4} - trial4.bias(4);
t4A1 = trial4.A{1};
t4A2 = trial4.A{2};
t4A3 = trial4.A{3};
t4A4 = trial4.A{4};

t7Z1 = trial7.Z{1} - trial7.bias(1);
t7Z2 = trial7.Z{2} - trial7.bias(2);
t7Z3 = trial7.Z{3} - trial7.bias(3);
t7Z4 = trial7.Z{4} - trial7.bias(4);
t7A1 = trial7.A{1};
t7A2 = trial7.A{2};
t7A3 = trial7.A{3};
t7A4 = trial7.A{4};

figure(10); clf(10);
hold on; grid on;
scatter(t4A1,t4Z1, '.')
scatter(t4A2,t4Z2, '.')
scatter(t4A3,t4Z3, '.')
scatter(t4A4,t4Z4, '.')

scatter(t7A1,t7Z1, '.')
scatter(t7A2,t7Z2, '.')
scatter(t7A3,t7Z3, '.')
scatter(t7A4,t7Z4, '.')
legend('t4_1','t4_2','t4_3','t4_4', 't7_1','t7_2','t7_3','t7_4')

A47 = [trial4.A_all, trial7.A_all];
Z47 = [trial4.Z_all, trial7.Z_all];

A6 = trial6.A_all;
Z6 = trial6.Z_all;

figure(11); clf(11);
scatter(A6,Z6)

%% Nonlinear Curve Fitting

A467 = [trial4.A_all, trial6.A_all, trial7.A_all]';
Z467 = [trial4.Z_all, trial6.Z_all, trial7.Z_all]';

% give more weight to smaller areas
% W = 2./A47;
W = Z467/max(Z467)+.02;

% Fit: 'Zsensing-Phantom'.
[xData, yData, weights] = prepareCurveData( A467, Z467, W );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [284.457039538282 -0.830703793768868];
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%% Plot fit with data

figure
subplot_er(2,1,1)
title('Model Validation Data','FontSize',13)
hold on
h = plot(xData, yData, '.b');
plot(fitresult,'g')
annotation('textbox',[.47 .8 .01 .01],'String','$Z = 347.8A^{-1.035}$','FitBoxToText','on', 'FontSize',12, 'BackgroundColor','w','Interpreter','latex');

% Label axes
% xlabel( 'mm^2' ,'FontSize',11);
ylabel( '\Omega','FontSize',11 );
grid on
ylim([0 3100])


% Plot residuals.
subplot_er( 2, 1, 2 );
title('Residuals','FontSize',13)
h = plot( fitresult, xData, yData, 'residuals');
% legend( h, 'Zsensing-Phantom - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'mm^2','FontSize',11);
ylabel( '\Delta\Omega','FontSize',11);
grid on


%% Validation
% use Trial 6

% figure; hold on
% for ii=1:4
%     c3 = trial6.bias(ii);
%     Aest(:,ii) = ( (trial6.Z_raw(:,ii)-c3) / c1 ) .^ (1/c2);
% 
%     scatter(Aest, trial6.Z_raw(:,ii)-c3)
% 
% end

Aest = ( (trial6.Z_all) / c1 ) .^ (1/c2);


figure;
subplot(1,2,1)
hold on;
scatter(trial6.A_all, Aest)
plot([0 1],[0 1],'r')
axis([0 3 0 3])

subplot(1,2,2)
hold on
scatter(Aest, trial6.Z_all)
xlim([0 3])

%%
ransac_max_distance = 20;
ransac_iter = 1;
for ii = 1:4
    ZsensingFitCurveRansac(trial4.Z{ii}', trial4.A{ii}', ransac_max_distance, ransac_iter, true);
end
% [idx_inliers, coeffs, coeffs_lin] = ZsensingPlotFitCurvesRansac(trial4.Z, trial4.A, ransac_max_distance, ransac_iter);

%%
%%%%%%%%%%%%%%%%%%%%
%%% Cadaver Data %%%
%%%%%%%%%%%%%%%%%%%%

directory_base = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data';

% Impedance data
load(fullfile(directory_base, 'Zsense_cadaver_data.mat'));

% LSTM network
load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data\Zsense_LSTMclassify-min.mat')

% Area pts from segmented CTs
voxel_size = 0.2; % [mm]
for i_trial = 2:5
    for i_ch = 1:4
        Apts{i_trial,i_ch} = voxel_size * importPts3D(fullfile(directory_base, 'TBone4', sprintf('Trial%i',i_trial), sprintf('Trial%i_Ch%i.txt', i_trial,i_ch)));
    end
end


%%
net = netZsensing;

cad(1).Z_raw = Z1_initialentry;
cad(2).Z_raw = Z2_CT1;
cad(3).Z_raw = Z3_CT2;
cad(4).Z_raw = Z4_CT3;
cad(5).Z_raw = Z5_CT4;


%% smoothing

span = 6; 

for ii = 1:length(cad)
    for jj=1:4
        cad(ii).Z(:,jj) = smoothdata(cad(ii).Z_raw(:,jj), 'movmean',span);
    end
end

%% Manually determine bias offsets for each channel
% cad(2).Z_bias = [];

%% LSTM Prediction

% find current minimums
for ii=1:length(cad)
    for i_ch=1:4
        for i_el=1:numel(cad(ii).Z_raw(:,i_ch))
            cad(ii).Z_min(i_el,i_ch) = min(cad(ii).Z_raw(1:i_el,i_ch));
        end
    end
end

% compute predictions
for ii = 1:length(cad)
    for i_ch = 1:4
        net.resetState;
        cad(ii).pred(:,i_ch) = classify(net, (cad(ii).Z_raw(:,i_ch) - cad(ii).Z_min(:,i_ch))');
    end
end

%% Final Z
% i.e. impedance when CT was acquired

n_avg = 500; % average the last N samples

for ii=1:5
    cad(ii).Z_final = mean(cad(ii).Z_raw(end-n_avg:end, :));
    fprintf('Z_final_CT%i = [%4.0f, %4.0f, %4.0f, %4.0f] \n', ii-1, cad(ii).Z_final);
end


%% Plot Z

figure
for i_trial = 1:5
    subplot_er(2,3,i_trial);
    title(sprintf('Trial %i', i_trial))
    plot(cad(i_trial).Z)
    ylabel('\Omega')
    ylim([2000 4500])
    legend
end

%% Compute Areas (segmented from CIP/IMPROVISE)

colors = distinguishable_colors(4);

figure
for i_trial = 2:5
    subplot_er(2,2,i_trial-1);
    hold on
    for i_ch = 1:4
        plane{i_trial,i_ch} = fitPlane(Apts{i_trial,i_ch});
        planar_pts{i_trial,i_ch} = projPointOnPlane(Apts{i_trial,i_ch}, plane{i_trial,i_ch});
        cad(i_trial).A(i_ch) = polygonArea3d(planar_pts{i_trial,i_ch});
        cad(i_trial).A2(i_ch) = polygonArea3d(Apts{i_trial,i_ch});
        drawPoint3d(Apts{i_trial,i_ch}, 'Marker','.', 'Color',colors(i_ch,:))
        drawPoint3d(planar_pts{i_trial,i_ch}, 'Marker','o', 'Color',colors(i_ch,:))
        h=fillPolygon3d(planar_pts{i_trial,i_ch}, colors(i_ch,:));
        axis vis3d
    end
end


%% Z vs A
figure; hold on
for i_trial = 2:5
    scatter(cad(i_trial).A2, cad(i_trial).Z_final)
end
legend('Location','sw')

%% 

