function [idx_inliers, coeffs, coeffs_lin, S, S_lin] = ZsensingFitCurveRansac(R, A, max_distance, iter, plot_flag)
%%

fitLineFcn  = @(X) polyfit(X(:,1),X(:,2),1); % fit function using polyfit
evalLineFcn = @(model, X) sum((X(:, 2) - polyval(model, X(:,1))).^2,2); % distance evaluation function


if nargin < 5
    plot_flag = false;

    if nargin < 4
        iter = 1;
        
        if nargin < 3
            max_distance = 0.5;
        end
    end
end


%% Linearize Data

% Put R = c1*A^c2 + c3 into linearized form (y=m*x+b)  ==>  ln(R-c3) = c2*ln(A) + ln(c1)
% y = ln(R-c3)
% m = c2
% x = ln(A)
% b = ln(c1)

c3_lin = min(R)-0.01; % [ohms] use min so R-c3 is always positive
Alin = log(A);
Rlin = log(R-c3_lin);
data = [Alin Rlin];

[model_lin, S_lin] = polyfit(Alin, Rlin, 1);
coeffs_lin = [exp(model_lin(2)), model_lin(1), c3_lin];


%% Run RANSAC to Determine Outliers

sample_size = 2;

data_new = data;
c3_new = c3_lin;

for ii=1:iter

    [~, idx_inliers_new{ii}] = ransac(data_new, fitLineFcn, evalLineFcn, sample_size, max_distance, 'MaxNumTrials',5e6, 'Confidence', 99.99);
    
    % refit line to inliers
    data_inliers = data_new(idx_inliers_new{ii},:);
    [model_inliers, S] = polyfit(data_inliers(:,1), data_inliers(:,2), 1); 
    c1(ii) = exp(model_inliers(2));
    c2(ii) = model_inliers(1);
    c3(ii) = c3_new;

    % remove outliers and re-linearize using new minimum R for c3
    idx_inliers_all = 1:length(R);
    for jj=1:ii
        idx_inliers_all = idx_inliers_all(idx_inliers_new{jj});
    end
    idx_inliers = logical(zeros(size(R)));
    idx_inliers(idx_inliers_all) = true;

    R_new = R(idx_inliers);
    A_new = A(idx_inliers);

    c3_new = min(R_new) - 0.01;
    data_new = [log(A_new), log(R_new-c3_new)];
%     [model_inliers, S] = polyfit(data_new(:,1), data_new(:,2), 1); % refit line to inliers



end

coeffs = [c1(end) c2(end) c3(end)];


%% Plot

if plot_flag
    figure; 
    
% Linearized
    subplot_er(2,2,1);
    hold on

    % Create scatter plot of linearized raw data
    scatter(Alin, Rlin);

    % mark outliers
    scatter(Alin(~idx_inliers), Rlin(~idx_inliers),'xr');

    % least-squares linear model fit
    Afit_lin = [min(Alin), max(Alin)];
    Rfit_lin = model_lin(1)*Afit_lin + model_lin(2);
    plot(Afit_lin, Rfit_lin, 'm-', 'LineWidth', 2) 

    ylim([0,10])
    legend('Raw Data', 'Outliers', 'Least squares fit');
    grid on
%     legend('Raw Data', 'Least squares fit');


% Robust Linearized
    subplot_er(2,2,3);
    hold on

    % Create scatter plot of linearized raw data
    scatter(data_inliers(:,1), data_inliers(:,2));

    Afit = [min(data_inliers(:,1)) max(data_inliers(:,1))];
    Rfit = model_inliers(1)*Afit + model_inliers(2);
    plot(Afit, Rfit, 'g-', 'LineWidth', 2)

    ylim([0,10])
    legend('Raw Data','Robust fit');
    grid on

%     subplot_er(1,2,1);
%     hold on
% 
%     % Create scatter plot of linearized raw data
%     scatter(Alin, Rlin);
% 
%     % mark outliers
%     scatter(Alin(~idx_inliers), Rlin(~idx_inliers),'xr');
% 
%     % least-squares linear model fit
%     Afit_lin = [min(Alin), max(Alin)];
%     Rfit_lin = model_lin(1)*Afit_lin + model_lin(2);
%     plot(Afit_lin, Rfit_lin, 'm-', 'LineWidth', 2) 
% 
%     % robust model fit
%     Afit = [min(data_inliers(:,1)) max(data_inliers(:,1))];
%     Rfit = model_inliers(1)*Afit + model_inliers(2);
%     plot(Afit, Rfit, 'g-', 'LineWidth', 2)
% 
%     ylim([0,10])
%     legend('Raw Data','Outliers','Least squares fit','Robust fit');


% Resistance vs Area
    subplot_er(2,2,[2,4]);
    hold on

    % transform from linearized form back to original: R = c1*A^c2 + c3;
    Rfit_linear = coeffs_lin(1)*A.^coeffs_lin(2) + coeffs_lin(3);
    Rfit_robust = coeffs(1)*A.^coeffs(2) + coeffs(3);

    % Create scatter plot of raw data
    scatter(A,R,'bo')

    % mark outliers
    scatter(A(~idx_inliers), R(~idx_inliers),'xr')

    % model fits
    [~,i_Asort] = sort(A);
    plot(A(i_Asort), Rfit_linear(i_Asort), 'm', 'LineWidth', 2)
    plot(A(i_Asort), Rfit_robust(i_Asort), 'g', 'LineWidth', 2)

    txt_lin    = ['R_{linear} = ', num2str(coeffs_lin(1), '%.1f'), '*A^{', num2str(coeffs_lin(2),3), '} + ', num2str(coeffs_lin(3),4)];
    txt_robust = ['R_{robust} = ', num2str(coeffs(1),     '%.1f'), '*A^{', num2str(coeffs(2),    3), '} + ', num2str(coeffs(3),    4)];
    txt_x = 1;
    txt_y = 0.5*(max(R)-min(R))+min(R);
    text(txt_x, txt_y,     txt_lin,    'FontSize', 14);
    text(txt_x, txt_y-200, txt_robust, 'FontSize', 14);

    legend('Raw Data','Outliers','Least Squares Fit','Robust Fit')
    xlabel('Area (mm^2)')
    ylabel('Resistance (\Omega)')

    xlim([0 max(A)])
    ylim([floor(min(R)/100)*100 max(R)])
    grid on
end

end