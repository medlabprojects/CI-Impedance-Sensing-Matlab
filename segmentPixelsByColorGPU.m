%% segmentPixelsByColor
%
%   rgb_image => [NxMx3] color image array
%   rgb_target_color => [1x3] desired color to segment, values must by 0->1 (e.g. red = [1 0 0])
%   hsv_tolerance => 'percentage' differences allowable (0->1); default = [0.05, 0.1, 0.1]
%
%   matching_pixels => [NxM] logical array of segmented pixels
%
%
%   Trevor Bruns
%   September 2019

function matching_pixels = segmentPixelsByColorGPU(rgb_image, rgb_target_color, hsv_tolerance)
    
    %% convert from RGB to HSV
    hsv_color = rgb2hsv(rgb_target_color);    % rgb2hsv seems to run much faster on CPU
    hsv_image = gpuArray(rgb2hsv(rgb_image)); % load image into GPU
    
    %% tolerance parameters
    if nargin < 3
    hsv_tolerance = [0.1, 0.1, 0.1]; % segment pixels with a hue +- this tolerance
    end
    
    hue_low  = hsv_color(1) - hsv_tolerance(1);
    hue_high = hsv_color(1) + hsv_tolerance(1);
    sat_low  = hsv_color(2) - hsv_tolerance(2);
    sat_high = hsv_color(2) + hsv_tolerance(2);
    val_low  = hsv_color(3) - hsv_tolerance(3);
    val_high = hsv_color(3) + hsv_tolerance(3);

    % check for hue wrap-around issues (bounds are 0->1)
    wrap_around = false;
    if hue_low < 0
        hue_high = hue_low + 1; % become high after wrap-around
        hue_low = hsv_color(1) + hsv_tolerance(1); % high becomes low
        wrap_around = true;
    elseif hue_high > 1
        hue_low = hue_high - 1; % become low after wrap-around
        hue_high = hsv_color(1) - hsv_tolerance(1); % low becomes high
        wrap_around = true;
    end
    
    % saturation and value are also bounded 0->1 but do not wrap-around
    if sat_low < 0
        sat_low = 0;
    end
    if sat_high > 1
        sat_high = 1;
    end
    
    if val_low < 0
        val_low = 0;
    end
    if val_high > 1
        val_high = 1;
    end
    
    %% segment matching pixels
    % hue
    if wrap_around
        % matches if between 0->low OR between high->1
        matching_pixels_gpu = (hsv_image(:,:,1) <= hue_low) | (hsv_image(:,:,1) >= hue_high);
    else
        % matches if between low->high
        matching_pixels_gpu = (hsv_image(:,:,1) >= hue_low) & (hsv_image(:,:,1) <= hue_high);
    end
    
    % saturation
    matching_pixels_gpu = matching_pixels_gpu & ((hsv_image(:,:,2) >= sat_low) & (hsv_image(:,:,2) <= sat_high));
    
    % value
    matching_pixels_gpu = matching_pixels_gpu & ((hsv_image(:,:,3) >= val_low) & (hsv_image(:,:,3) <= val_high));
            
    %% gather back from GPU
    matching_pixels = gather(matching_pixels_gpu);
end