%% Load data

load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data\Zsense_dissertation-data_trained-net.mat');

%% Prediction using LSTM network

bin_A_net = [];
for ii = 1:4
    bin_A_net = [bin_A_net, trial6.A_bin_pred{ii}];
end


%% Prediction using model fit
% Model:  
% (R-Rmin) = 347.8*A^(-1.035)
% ((R-Rmin)/347.8)^(-1/1.035) = A

c1 = 347.8;
c2 = -1.035;

% estimate for A
A_est = (trial6.Z_all/c1).^(1/c2);

% determine bin for each estimate
bin_A_model = discretize(A_est, [bin_edges(1:3) inf]);

% true bin
bin_A = discretize(trial6.A_all, [bin_edges(1:3) inf]);

% convert to bin names
bin_A_model = categorical(bin_A_model, 1:length(bin_names), bin_names); 
bin_A = categorical(bin_A, 1:length(bin_names), bin_names);


%% Plot: prediction vs time

figure; hold on
plot(bin_A, 'k--', 'LineWidth',1.5)
plot(bin_A_model, 'b', 'LineWidth',1.2)
plot(bin_A_net, 'm', 'LineWidth',1.2)
xlabel('Time (# samples)')
legend('Actual','Model Prediction', 'LSTM Prediction')


%% Confusion Matrix

figure;
subplot_er(1,2,1);
h_cm(1) = confusionchart(bin_A, bin_A_model,'RowSummary','row-normalized','ColumnSummary','column-normalized');
h_cm(1).Title = 'Nonlinear Model';
subplot_er(1,2,2);
h_cm(2) = confusionchart(bin_A, bin_A_net,'RowSummary','row-normalized','ColumnSummary','column-normalized');
h_cm(2).Title = 'LSTM Network';

%% Plot: true area vs time
% figure; hold on
% title('Predicted vs Measured Areas')
% for i_class = 1:numClasses        
%     % shade the region corresponding to this prediction
%     fill([0 0 length(YPred{i_ch}) length(YPred{i_ch})],...
%          [bin_edges(i_class) bin_edges(i_class+1) bin_edges(i_class+1) bin_edges(i_class)],...
%          class_colors(i_class,:), 'FaceAlpha',0.1, 'EdgeColor','none')
% 
%     % find all samples with this prediction
%     idx = find(YPred{i_ch} == bin_names(i_class));
% 
%     % plot all samples with this prediction
%     plot(idx,trial6.A{i_ch}(1,idx),'.', 'Color',class_colors(i_class,:), 'MarkerSize',10)
% end
%     
% ylim([0 2.5])
% ylabel('Area (mm^2)', 'FontSize',10)
% end
% linkaxes(ax_area,'x')
% xlabel("Time Step")
% hold off


