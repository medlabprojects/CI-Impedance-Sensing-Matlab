function [sync_score, Z_sync] = ZsensingCheckSync(Z_start, Z_end, areas, Z, varargin)
%%
%{

areas   => [Nx1] array of areas {mm^2} segmented from video
Z_raw   => [Mx1] array of all impedance measurements
Z_start => start index of Z_raw corresponding to start of areas
Z_end   => end index of Z_raw corresponding to end of areas

score   => arbitrary 'score' for how well the data matches, lower is better
Z_sync  => synced and interpolated Z values

%}

%% Parse Inputs

% input validation functions
isPosInt      = @(x) (x>=0) && (abs(round(x)-x) < eps); % verifies that x is a positive integer (works even if user doesn't specify since e.g. 'x=4' defaults to a double)
validDataSize = @(x) (size(x,2) == 1) && (size(x,1) > size(x,2)); % vertical array of samples (i.e. more rows than columns)
validMethods  = @(x) any(validatestring(x,{'linear','ransac'}));

% defaults
default_method = 'linear';
default_plot   = false;

% define parser
p = inputParser;
addRequired(p, 'Z_start', isPosInt);
addRequired(p, 'Z_end',   isPosInt);
addRequired(p, 'areas',   validDataSize);
addRequired(p, 'Z',       validDataSize);
addParameter(p, 'Plot',   default_plot, @islogical);
addParameter(p, 'Method', default_method, validMethods);

% parse and set vars
parse(p,Z_start,Z_end,areas,Z,varargin{:});
plot_flag = p.Results.Plot;
method    = p.Results.Method;


% rescale to match video frames
Z_sync_ind = linspace(Z_start, Z_end, length(areas))';               % rescale to get (fractional) indices corresponding to each area
Z_sync = interp1(Z_start:Z_end, Z(Z_start:Z_end), Z_sync_ind);  % interpolate at these indices to get impedances

% fit curve
if strcmp(method,'linear')
    [coeffs, S] = ZsensingFitCurve(Z_sync, areas, plot_flag); % linearized
elseif strcmp(method,'ransac')
    [~, coeffs, ~, S] = ZsensingFitCurveRansac(Z_sync, areas, 1, 1, plot_flag); % RANSAC
else
    error('invalid method specified')
end

sync_score = S.normr;


% compute norm of residuals
% Z_fit = coeffs(1)*areas.^coeffs(2) + coeffs(3);
% normr = sqrt(sum( (Z_fit-Z_sync).^2 ));
% sync_score = normr;



if plot_flag

    Z_fit = coeffs(1)*areas.^coeffs(2) + coeffs(3);


    % normalize to make comparison easier
%     areas_inv_sq = areas.^-1.2; % relationship is ~proportional to 1/A^-c
%     areas_norm = (areas_inv_sq-min(areas_inv_sq)) / max(areas_inv_sq-min(areas_inv_sq)); % normalize so range is [0,1]
    Z_fit_norm = (Z_fit-min(Z_fit)) / max(Z_fit-min(Z_fit));
    Z_norm = (Z_sync-min(Z_sync)) / max(Z_sync-min(Z_sync));
%     norm_diff = areas_norm(:) - Z_norm(:); 

    % 
    figure; 
%     subplot(2,1,1)
    hold on;
    plot(Z_sync,'b')
    plot(Z_fit, 'k')
%     plot(Z_norm,'b')
%     plot(Z_fit_norm, 'k')
%     plot(areas_norm, 'k')
%     plot(norm_diff, 'r')
    title(sprintf('Z Sync (''score'' = %.4f)', sync_score));
    legend('Z_{sync}', 'Z_{fit}', 'Location','nw')
    grid minor

% 	subplot(2,1,2)
%     scatter(xlin, ylin) % Create scatter plot of linearized data
%     hold on
%     xfit = [min(xlin), max(xlin)];
%     yfit = p(1)*xfit + p(2);
%     plot(xfit, yfit, 'm-', 'LineWidth',2)
end

end