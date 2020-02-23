%% Import Impedance Data
% need to create Nx1 array called Zraw

directory_name = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2020-01-20';

% load(fullfile(directory_name, 'Zsense_2020-01-20_phantom_Flex20_EA-2-3-5_trial1.mat'));

% trim calibration date (due to 'shift' of some kind) and compute mean
trim_ind = 750;
Z23_offset = mean(Z_2_3_cal(trim_ind:end));
Z35_offset = mean(Z_3_5_cal(trim_ind:end));

% apply offsets
Z23 = Z_2_3 - Z23_offset;
Z35 = Z_3_5 - Z35_offset;

% create Zraw
Z_raw = Z23;

%% Load Trial 2

directory_name = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2020-01-20';
load(fullfile(directory_name, 'Zsense_2020-01-20_phantom_Flex20_EA-2-3_trial2.mat'));
Z_raw = Z23_trial2; % no calibration
% Z_raw = Z23_trial2 - 2400; % something wrong with calibration data (saline evaporation from sitting out?)



%% smooth Z_raw
smooth_span = 30; % will depend on the acquisition rate
Z_smooth = smooth(Z_raw, smooth_span, 'rloess');
% Z_smooth = smooth(Z_raw, smooth_span, 'sgolay',3);

figure; hold on;
plot(Z_raw,'k')
plot(Z_smooth,'r')



%% Import Video

% filename_video = 'Zsense_2020-01-20_phantom_Flex20_EA-2-3-5_trial1_tracked';
filename_video = 'Zsense_2020-01-20_phantom_Flex20_EA-2-3_trial2_tracked';
vid = VideoReader( fullfile(directory_name, strcat(filename_video, '.mp4')) );

% image scale- determined by measuring width of phantom (10mm) in pixels
hf_scale = imshow(vid.readFrame, 'border','tight');
[scale.x, scale.y] = ginput(2);

% pixel_scale = 821/(10e-3); % [pixels/meter] for trial 1
pixel_scale = abs(scale.x(2)-scale.x(1))/(10e-3); % [pixels/meter] 

pixel_area = (1e3 / pixel_scale)^2; % [mm^2/pixel]



%% Segment Shaded Region of Interest (ROI) From Video

shaded_color_rgb = [1 0 0]; % red
hsv_tol = [0.03 0.1 0.1]; % allowable 'percentage' difference 

% trial1
% start_time = 9.25; % [s]
% start_time = 23; % [s]
% end_time   = 82;   % [s]

% trial2
start_time = 0; % [s]
end_time   = 22.4;   % [s]

start_frame = floor(start_time * vid.FrameRate);
end_frame = floor(end_time * vid.FrameRate);
frames = start_frame:end_frame;
total_frames = length(frames);

tic;
spmd
    % divide frames among individual workers
    my_frames = frames( discretize(frames, numlabs) == labindex );
    my_start_frame = my_frames(1);

    % go to first frame (NOTE: Assumes constant frame rate!)
    my_vid = vid; % each worker needs its own vidobj copy
    my_vid.CurrentTime = (1/my_vid.FrameRate) * my_start_frame; 

    % segment
    my_areas = zeros(length(my_frames),1);
    for ii = 1:length(my_frames)
        
        if ~my_vid.hasFrame % first ensure we still have frames to read
            break;
        end

        shaded_pixels = segmentPixelsByColorGPU(my_vid.readFrame, shaded_color_rgb);
        my_areas(ii) = sum(shaded_pixels(:)) * pixel_area;
    end
end

% assemble results into single variable
areas_raw = zeros(total_frames,1);
for ii = 1:length(my_areas)
    areas_raw(my_frames{ii}-start_frame+1) = my_areas{ii};
end
toc

% smooth
smooth_span_vid = 0.3 * vid.FrameRate; % smooth over 0.3 second span
areas = smooth(areas_raw, smooth_span_vid, 'sgolay', 4)';

figure; hold on;
plot(areas_raw,'k')
plot(areas,'r')



%% Sync Data/Video and Verify

% pick lower and upper bounds from plot
figure; plot(Z_smooth);
[bounds,~] = ginput(4);
lb = round([bounds(1), bounds(3)]); % must be integers
ub = round([bounds(2), bounds(4)]);

% lb = [2800, 9400];
% ub = [3500, 9900];

% use genetic algorithm to 'fine-tune' the fit (note: fmincon cant be used for integer values)
sync_fun = @(x)ZsensingCheckSync(x(1), x(2), areas, Z_smooth);
opts = optimoptions('ga','PlotFcn',@gaplotbestf);
x = ga(sync_fun,2,[],[],[],[],lb,ub,[],[1,2],opts)

% use indices found to compute Z_sync and create verification plot
[sync_score, Z_sync] = ZsensingCheckSync(x(1), x(2), areas, Z_smooth, true);



%% Setup Figure

makeVideo = false;
filename_new_video = strcat(filename_video, '_Matlab.mp4');

dpi_scale = 96/144; % matlab assumes 96dpi on Windows, but Dell M3800 is actually 144dpi
% figWidth = 1280*dpi_scale;
% figHeight = 720*dpi_scale;
figWidth = 1920*dpi_scale;
figHeight = 1080*dpi_scale;
gap = figHeight/100;

Z = Z_sync;

hFig = figure('units','pixels','Position',[100 10 figWidth figHeight],'MenuBar','none');

% video frame
hVid   = axes('units','pixels', 'Position',[0, 0, floor(vid.Width*figHeight/vid.Height), figHeight]);
hImage = imshow(vid.readFrame,'border','tight', 'Parent',hVid);

% Z vs time
hZ = axes('units','pixels');
hold(hZ, 'on')
Z_times = (1:length(Z))' / vid.FrameRate; % [s]
ZaxisLimits = [0, Z_times(end), min(Z)-50, max(Z)+100];
axis(hZ, ZaxisLimits);
hZ_line = animatedline(Z_times(1), Z(1));
hZ_marker = line(hZ, Z_times(1), Z(1), 'Color', 'm', 'LineStyle','none', 'Marker','o', 'MarkerSize',12);
xlabel(hZ, 'Time [s]');
ylabel(hZ, 'Z [\Omega]')

% Z vs Area
hArea  = axes('units','pixels');
hold(hArea, 'on')
AreaAxisLimits = [0, 1.05*max(areas), ZaxisLimits(3:4)];
axis(hArea, AreaAxisLimits);
hArea_line = animatedline(hArea, areas(1), Z(1), 'Color', 'r', 'LineStyle','none', 'Marker','x', 'MarkerSize',5);
hArea_marker = line(hArea, areas(1), Z(1), 'Color', 'g', 'LineStyle','none', 'Marker','.', 'MarkerSize',30);
ylabel(hArea, 'Z [\Omega]');
xlabel(hArea, 'Area between electrodes [mm^2]');

% position each plot
margins = max([hZ.TightInset; hArea.TightInset]); % pick largest
hArea.Position = [hVid.Position(3)+margins(1)+gap/2, margins(2)+gap,             figWidth-hVid.Position(3)-margins(1)-margins(3)-1.5*gap, figHeight/2-margins(2)-margins(4)-2*gap];
hZ.Position    = [hVid.Position(3)+margins(1)+gap/2, figHeight/2+margins(2)+gap, figWidth-hVid.Position(3)-margins(1)-margins(3)-1.5*gap, figHeight/2-margins(2)-margins(4)-2*gap];

% set up VideoWriter
if makeVideo
    videoFWriter = vision.VideoFileWriter(fullfile(directory_name, filename_new_video), 'FileFormat','MPEG4', 'FrameRate',vid.FrameRate, 'Quality',80);
end


%% plot/create video

vid.CurrentTime = start_time; % go to first frame
currentFrame = 0; % reset count

tic;
while (hasFrame(vid) && (currentFrame < total_frames))

    % read next frame
    currentFrame = currentFrame + 1;
    vidFrame = readFrame(vid);


    % Update Plots

    % current frame image
    set(hImage,'CData', vidFrame); % faster than calling imshow


    % Time vs Z plot
    delete(hZ_marker);
    addpoints(hZ_line,   Z_times(currentFrame), Z(currentFrame));
    hZ_marker = line(hZ, Z_times(currentFrame), Z(currentFrame), 'Color', 'm', 'LineStyle','none', 'Marker','o', 'MarkerSize',12);

    % Area vs Z plot
    delete(hArea_marker);
    addpoints(hArea_line,      areas(currentFrame), Z(currentFrame));
    hArea_marker = line(hArea, areas(currentFrame), Z(currentFrame), 'Color', 'g', 'LineStyle','none', 'Marker','.', 'MarkerSize',30);


    if makeVideo
        drawnow
        frame = getframe(hFig);
        videoFWriter(frame.cdata);
    else
        pause(1/vid.FrameRate);
    end

end

if makeVideo
    release(videoFWriter);
end
toc