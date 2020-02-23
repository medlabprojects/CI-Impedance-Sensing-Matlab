function [coeffs, S] = ZsensingFitCurve(R, A, plot_flag)
%%
%   R = [Nx1]
%   A = [Nx1]
%   
%   coeffs = [c1 c2 c3]


if nargin < 3
    plot_flag = false;
end

%% Linearize Data

% Put R = c1*A^c2 + c3 into linearized form (y=m*x+b)  ==>  ln(R-c3) = c2*ln(A) + ln(c1)
% y = ln(R-c3)
% m = c2
% x = ln(A)
% b = ln(c1)

c3 = min(R)-0.01; % [ohms] use min so R-c3 is always positive
Alin = log(A);
Rlin = log(R-c3);

[model, S] = polyfit(Alin, Rlin, 1);
coeffs = [exp(model(2)), model(1), c3];


%% Plot

if plot_flag
    figure; 
    
% Linearized
    subplot_er(1,2,1);
    hold on

    % Create scatter plot of linearized raw data
    scatter(Alin, Rlin);

    % least-squares linear model fit
    Afit_lin = [min(Alin), max(Alin)];
    Rfit_lin = model(1)*Afit_lin + model(2);
    plot(Afit_lin, Rfit_lin, 'm-', 'LineWidth', 2) 

%     ylim([0,10])
    legend('Raw Data', 'Least squares fit');


% Resistance vs Area
    subplot_er(1,2,2);
    hold on

    % transform from linearized form back to original: R = c1*A^c2 + c3;
    Rfit = coeffs(1)*A.^coeffs(2) + coeffs(3);

    % Create scatter plot of raw data
    scatter(A,R,'bo')

    % model
    [~,i_Asort] = sort(A);
    plot(A(i_Asort), Rfit(i_Asort), 'm', 'LineWidth', 2)

    txt_lin    = ['R_{linear} = ', num2str(coeffs(1), '%.1f'), '*A^{', num2str(coeffs(2),3), '} + ', num2str(coeffs(3),4)];
    txt_x = 1;
    txt_y = 0.5*(max(R)-min(R))+min(R);
    text(txt_x, txt_y,     txt_lin,    'FontSize', 14);

    legend('Raw Data','Least Squares Fit')
    xlabel('Area (mm^2)')
    ylabel('Resistance (\Omega)')

    xlim([0 max(A)])
    ylim([floor(min(R)/100)*100 max(R)])
end

end