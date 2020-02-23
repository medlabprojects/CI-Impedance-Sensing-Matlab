function [idx_Z_sync, Z_sync, sync_scores] = ZsensingSyncData(areas, Z, varargin)
%% Uses genetic algorithm to determine the indices of Z that best sync with start/end of areas
%
%   areas => [NxM] N rows of measurements of M individual areas (e.g. 200 measurements using 4 channels => [200x4]
%   Z     => [NxM] impedance measurements corresponding to each area
%   
%   Optional:
%   'Channel' => [scalar] which channel to use when syncing (default = M, the last channel)
%   'Bounds'  => [1x4]    indice ranges in Z to use for syncing [lb_min, lb_max, ub_min, ub_max]
%   'Plot'    => [bool]   generate verification plots


%% Parse Inputs

% input validation functions
isPosInt      = @(x) (x>=0) && (abs(round(x)-x) < eps); % verifies that x is a positive integer (works even if user doesn't specify since e.g. 'x=4' defaults to a double)
validDataSize = @(x) (length(size(x)) == 2) && (size(x,1) > size(x,2)); % 2D matrix with each row being a sample (i.e. more rows than columns)
validZSize    = @(x) size(x,2)==size(areas,2); % must have the same number of channels/areas per sample
validChannel  = @(x) isPosInt(x) && (x <= size(areas,2));
validBounds   = @(x) (size(x,2)==size(Z,2)) && (0<x(1)) && (x(1)<x(2)) && (x(2)<x(3)) && (x(3)<x(4)) && (x(4)<=size(Z,1));
validMethods  = @(x) any(validatestring(x,{'brute','ga'}));

% defaults
default_channel = size(Z,2); % use last channel (most likely to have largest spread of values)
default_bounds = [0 0 0 0];  % not actually used, will prompt user with ginput
default_method = 'ga';

% define parser
p = inputParser;
addRequired(p, 'areas', validDataSize);
addRequired(p, 'Z', validZSize);
addParameter(p, 'Channel', default_channel, validChannel);
addParameter(p, 'Bounds', default_bounds, validBounds);
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'Method', default_method, validMethods);

% parse and set vars
parse(p,areas,Z,varargin{:});
i_ch      = p.Results.Channel;
plot_flag = p.Results.Plot;
method    = p.Results.Method;

if sum(p.Results.Bounds) < 1 % i.e. default of all zeros
    % have user select lower and upper bounds from plot

    hf_areas = figure; 
    plot(areas(:,i_ch),'k') % for comparison

    hf_Z = figure; 
    plot(Z(:,i_ch));
    [bounds,~] = ginput(4);
    lb = round([bounds(1), bounds(3)]); % must be integers
    ub = round([bounds(2), bounds(4)]);
    close(hf_areas)
    close(hf_Z)
else
    lb = [p.Results.Bounds(1), p.Results.Bounds(3)];
    ub = [p.Results.Bounds(2), p.Results.Bounds(4)];
end


%% Sync Data (Genetic Algorithm)

% use genetic algorithm to 'fine-tune' the fit (note: fmincon cant be used for integer values)
if strcmp(method,'ga')
    sync_fun = @(x)ZsensingCheckSync(x(1), x(2), areas(:,i_ch), Z(:,i_ch), 'Method','linear');
    opts = optimoptions('ga','PlotFcn',@gaplotbestf);
    idx_Z_sync = ga(sync_fun,2,[],[],[],[],lb,ub,[],[1,2],opts);
end


%% Sync Data (Brute Force)

if strcmp(method,'brute')

    % start with a course search
    span = 3;
    guess_low = lb(1):span:ub(1);
    guess_high = lb(2):span:ub(2);

    guess = zeros(length(guess_low)*length(guess_high), 2);
    count = 0;
    for ii=1:length(guess_low)
        for jj=1:length(guess_high)
            count = count + 1;
            guess(count,:) = [guess_low(ii), guess_high(jj)];
        end
    end

    A = areas(:,i_ch);
    R = Z(:,i_ch);
    normr = zeros(length(guess),1);
    parfor ii=1:length(guess)
        normr(ii) = ZsensingCheckSync(guess(ii,1), guess(ii,2), A, R, 'Method','linear');
    end

    [~,idx_best_guess] = min(normr);
    best_guess = guess(idx_best_guess,:);


    % fine search around the best guess so far
    range = 30;
    guess_low = (best_guess(1)-range/2):(best_guess(1)+range/2);
    guess_high = (best_guess(2)-range/2):(best_guess(2)+range/2);

    if guess_low(1) < lb(1)
        guess_low = lb(1):(best_guess(1)+range/2);
    end
    if guess_low(end) > ub(1)
        guess_low = (best_guess(1)-range/2):ub(1);
    end
    if guess_high(1) < lb(2)
        guess_high = lb(2):(best_guess(2)+range/2);
    end
    if guess_high(end) > ub(2)
        guess_high = (best_guess(2)-range/2):ub(2);
    end

    guess = zeros(length(guess_low)*length(guess_high), 2);
    count = 0;
    for ii=1:length(guess_low)
        for jj=1:length(guess_high)
            count = count + 1;
            guess(count,:) = [guess_low(ii), guess_high(jj)];
        end
    end

    normr = zeros(length(guess),1);
    parfor ii=1:length(guess)
        normr(ii) = ZsensingCheckSync(guess(ii,1), guess(ii,2), A, R, 'Method','linear');
    end

    [~,idx_best_guess] = min(normr);
    idx_Z_sync = guess(idx_best_guess,:);
end


%% Interpolate
% use indices found to interpolate Z corresponding to each area and (optionally) create verification plots
Z_sync = zeros(size(areas));
for ii = 1:size(Z,2)
    [sync_scores(ii), Z_sync(:,ii)] = ZsensingCheckSync(idx_Z_sync(1), idx_Z_sync(2), areas(:,ii), Z(:,ii), 'Plot',plot_flag, 'Method','ransac');
end

end