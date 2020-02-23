function [idx_inliers, coeffs, coeffs_lin] = ZsensingPlotFitCurvesRansac(Z, A, ransac_max_distance, ransac_iter)
%{
Designed for 4-channel data!!
Z = [Nx4]
A = [Nx4]

Fits curves to data and creates a figure with 2x2 subplots 


Trevor Bruns
February 2020

%}


%%

n_channels = 4;

Z_min = min(Z(:));
Z_max = max(Z(:));

figure;

idx_inliers = cell(1,n_channels);
coeffs = cell(1,n_channels);
coeffs_lin = cell(1,n_channels);

for ii = 1:n_channels
    subplot_er(2,2,ii);    
    hold on

    [idx_inliers{ii}, coeffs{ii}, coeffs_lin{ii}] = ZsensingFitCurveRansac(Z(:,ii), A(:,ii), ransac_max_distance, ransac_iter, false);

    % R = c1*A^c2 + c3
    Rfit_linear = coeffs_lin{ii}(1) * A(:,ii) .^ coeffs_lin{ii}(2) + coeffs_lin{ii}(3);
    Rfit_robust = coeffs{ii}(1)     * A(:,ii) .^ coeffs{ii}(2)     + coeffs{ii}(3);

    scatter(A(:,ii), Z(:,ii), 'bo') % Create scatter plot of raw data
    scatter(A(~idx_inliers{ii},ii), Z(~idx_inliers{ii},ii),'xr') % mark outliers

    [~,i_Asort] = sort(A(:,ii));
    plot(A(i_Asort,ii), Rfit_linear(i_Asort), 'm', 'LineWidth', 2)
    plot(A(i_Asort,ii), Rfit_robust(i_Asort), 'g', 'LineWidth', 2)

    txt_lin    = ['R_{linear} = ', num2str(coeffs_lin{ii}(1), '%.1f'), '*A^{', num2str(coeffs_lin{ii}(2),3), '} + ', num2str(coeffs_lin{ii}(3),4)];
    txt_robust = ['R_{robust} = ', num2str(coeffs{ii}(1),     '%.1f'), '*A^{', num2str(coeffs{ii}(2),    3), '} + ', num2str(coeffs{ii}(3),    4)];
    txt_x = 1;
    txt_y = 0.5*(Z_max-Z_min)+Z_min;
    text(txt_x, txt_y,     txt_lin,    'FontSize', 14);
    text(txt_x, txt_y-200, txt_robust, 'FontSize', 14);

    legend('Raw Data','Outliers','Least Squares Fit','Robust Fit')
    xlabel('Area (mm^2)')
    ylabel('Resistance (\Omega)')
    grid on

    xlim([0 max(A(:))])
    ylim([floor(Z_min/100)*100 Z_max])
end

end