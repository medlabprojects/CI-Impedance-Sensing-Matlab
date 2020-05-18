function [areas, polygons] = ZsensingSegmentVideoParallel(vid, start_time, end_time, pixels_per_meter)
%% 
%   - Segment processed (via After Effects) video
%   - Create Regions Between Electrodes/Modiolus
%   - Compute Areas
%
%   Trevor Bruns
%   February 2020
%
%
%%

if nargin < 4
    % image scale- determined by measuring width of phantom (10mm) in pixels
    hf_scale = figure;
    imshow(vid.readFrame, 'border','tight');
    hold on
    text(60,850, 'Calibration- Select a point on the left and right edge of phantom', 'Color','m', 'FontSize',14)
    [scale.x, scale.y] = ginput(2);

    pixel_scale = abs(scale.x(2)-scale.x(1))/(10e-3); % [pixels/meter] 
    pixel_area = (1e3 / pixel_scale)^2;               % [mm^2/pixel]

    close(hf_scale)

    if nargin < 3
        end_time = vid.Duration;
        
        if nargin < 2
            start_time = 0.01;
        end
    end

else
    pixel_area = (1e3/pixels_per_meter)^2; % [mm^2/pixel]
end

start_frame = floor(start_time * vid.FrameRate);
end_frame = floor(end_time * vid.FrameRate);
frames = start_frame:end_frame;
total_frames = length(frames);


%% Start parallel pool

if isempty(gcp('nocreate'))
    parpool('local', 4, 'IdleTimeout',Inf);
end


%% Segmentation

spmd
    % divide frames among individual workers
    n = floor(total_frames/numlabs);
    if labindex < numlabs
        my_frames = frames( (n*(labindex-1)+1) : (n*labindex) );
    else
        my_frames = frames( (n*(labindex-1)+1) : end ); % last worker takes any remainder
    end
    my_start_frame = my_frames(1);

    % go to first frame (NOTE: Assumes constant frame rate!)
    my_vid = vid; % each worker needs its own vidobj copy
    my_vid.CurrentTime = (1/my_vid.FrameRate) * my_start_frame; 

    my_regions(length(my_frames)) = struct();
    for ii = 1:length(my_frames)
        
        if ~my_vid.hasFrame % first ensure we still have frames to read
            break;
        end

        % Segment current frame to find centers of electrode markers and outlines of modiolus/EA masks
        [my_regions(ii).polygons, my_regions(ii).areas] = ZsensingFindRegions(my_vid.readFrame);
    end
end

% assemble worker vars into single matrices
regions = struct('polygons',{},'areas',{});
for ii = 1:length(my_regions)
    regions = [regions, my_regions{ii}];    
end

% outputs
areas = abs([regions.areas]') * pixel_area;
polygons = {regions.polygons}';

end