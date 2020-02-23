%% Import Impedance Data

directory_name = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2020-01-23';

% load(fullfile(directory_name, 'Zsense_2020-01-24_phantom_Flex24_5x20_EA-1-5_trial4.mat'));
% Z_raw = Z_trial4;

load(fullfile(directory_name, 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial7.mat'));
Z_raw = Z_trial7(1500:5500, :); % approximate range of trimmed video

% load(fullfile(directory_name, 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6.mat'));
% Z_raw = Z_trial6(1400:end, :); % approximate range of trimmed video

n_channels = size(Z_raw,2);
colors = distinguishable_colors(n_channels);


%% Smooth Z_raw

Z = Z_raw;
% smooth_span = 20; % will depend on the acquisition rate
% 
% Z = zeros(size(Z_raw));
% for ii=1:n_channels
% %     Z(:,ii) = smooth(Z_raw(:,ii), smooth_span, 'rloess');
%     Z(:,ii) = smooth(Z_raw(:,ii), smooth_span, 'sgolay',3);
% end
% 
% figure; hold on;
% plot(Z)
% plot(Z_raw,'.')


%% Import Videos

% filename_video     = 'Zsense_2020-01-24_phantom_Flex24_5x20_EA-1-5_trial4_crop-overlay'; % video with tracked overlays
% filename_video_raw = 'Zsense_2020-01-24_phantom_Flex24_5x20_EA-1-5_trial4_crop';         % same as ^ but without overlays
filename_video     = 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial7_crop-overlay2'; % video with tracked overlays
filename_video_raw = 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial7_crop';         % same as ^ but without overlays
% filename_video     = 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6_crop-overlay'; % video with tracked overlays
% filename_video_raw = 'Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6_crop';         % same as ^ but without overlays

vid     = VideoReader( fullfile(directory_name, strcat(filename_video, '.mp4')) );
vid_raw = VideoReader( fullfile(directory_name, strcat(filename_video_raw, '.mp4')) );


%% image scale- determined by measuring width of phantom (10mm) in pixels
hf_scale = figure;
imshow(vid.readFrame, 'border','tight');
hold on
text(60,850, 'Calibration- Select a point on the left and right edge of phantom', 'Color','m', 'FontSize',14)
[scale.x, scale.y] = ginput(2);

pixels_per_meter = abs(scale.x(2)-scale.x(1))/(10e-3); % [pixels/meter] 
pixel_area = (1e3 / pixels_per_meter)^2;               % [mm^2/pixel]

close(hf_scale)


%% Create Regions Between Electrodes/Modiolus and Compute Areas

start_time = 0.01;  % [s]
end_time   = vid.Duration; % [s]
% start_time = 30;  % [s]
% end_time   = 35; % [s]

tic;
[areas_raw, polygons] = ZsensingSegmentVideoParallel(vid, start_time, end_time, pixels_per_meter);
toc


%% Smooth Areas

smooth_span_vid = round(0.3 * vid.FrameRate); % smooth over X second span
areas = zeros(size(areas_raw));

figure; hold on;
for ii = 1:size(areas_raw,2)
    areas(:,ii) = smooth(areas_raw(:,ii), smooth_span_vid, 'sgolay', 4)';
    h_areas_raw(ii) = plot(areas_raw(:,ii), '.', 'Color',[colors(ii,:) 0.5]);
    plot(areas(:,ii), 'Color',colors(ii,:));
end
legend(h_areas_raw,'1','2','3','4')


%% Sync Data/Video and Verify

% [idx_Z_sync, Z_sync, sync_scores] = ZsensingSyncData(areas, Z, 'Channel',4, 'Method','ga', 'Plot',true);
[idx_Z_sync, Z_sync, sync_scores] = ZsensingSyncData(areas, Z, 'Channel',4, 'Method','brute', 'Bounds', [470 570 3750 3850], 'Plot',true); % trial 7
% [idx_Z_sync, Z_sync, sync_scores] = ZsensingSyncData(areas, Z, 'Channel',4, 'Method','brute', 'Bounds', [400 500 3325 3425], 'Plot',true); % trial 6
% [idx_Z_sync, Z_sync, sync_scores] = ZsensingSyncData(areas, Z, 'Channel',3, 'Method','brute', 'Bounds', [575 675 6575 6675], 'Plot',true); % trial 4

Z_min = min(Z_sync(:));
Z_max = max(Z_sync(:));


%% Compute LSTM Model Predictions for each channel
% load("Zsense_LSTMclassify.mat") % load 'net'
% for ii=1:4
%     pred(ii,:) = classify(net, trial6.Z_bias{ii});
% end


%% Plot Fit Curves

ransac_max_distance = 1.2;
ransac_iter = 2;
[idx_inliers, coeffs, coeffs_lin] = ZsensingPlotFitCurvesRansac(Z_sync, areas, ransac_max_distance, ransac_iter);


% %% Plot Fit Curves
% 
% ransac_max_distance = 1.2;
% ransac_iter = 2;
% 
% figure;
% 
% idx_inliers = cell(1,n_channels);
% coeffs = cell(1,n_channels);
% coeffs_lin = cell(1,n_channels);
% 
% for ii = 1:n_channels
%     subplot_er(2,2,ii);    
%     hold on
% 
%     [idx_inliers{ii}, coeffs{ii}, coeffs_lin{ii}] = ZsensingFitCurveRansac(Z_sync(:,ii), areas(:,ii), ransac_max_distance, ransac_iter, false);
% 
%     % R = c1*A^c2 + c3
%     Rfit_linear = coeffs_lin{ii}(1) * areas(:,ii) .^ coeffs_lin{ii}(2) + coeffs_lin{ii}(3);
%     Rfit_robust = coeffs{ii}(1)     * areas(:,ii) .^ coeffs{ii}(2)     + coeffs{ii}(3);
% 
%     scatter(areas(:,ii), Z_sync(:,ii), 'bo') % Create scatter plot of raw data
%     scatter(areas(~idx_inliers{ii},ii), Z_sync(~idx_inliers{ii},ii),'xr') % mark outliers
% 
%     [~,i_Asort] = sort(areas(:,ii));
%     plot(areas(i_Asort,ii), Rfit_linear(i_Asort), 'm', 'LineWidth', 2)
%     plot(areas(i_Asort,ii), Rfit_robust(i_Asort), 'g', 'LineWidth', 2)
% 
%     txt_lin    = ['R_{linear} = ', num2str(coeffs_lin{ii}(1), '%.1f'), '*A^{', num2str(coeffs_lin{ii}(2),3), '} + ', num2str(coeffs_lin{ii}(3),4)];
%     txt_robust = ['R_{robust} = ', num2str(coeffs{ii}(1),     '%.1f'), '*A^{', num2str(coeffs{ii}(2),    3), '} + ', num2str(coeffs{ii}(3),    4)];
%     txt_x = 1;
%     txt_y = 0.5*(Z_max-Z_min)+Z_min;
%     text(txt_x, txt_y,     txt_lin,    'FontSize', 14);
%     text(txt_x, txt_y-200, txt_robust, 'FontSize', 14);
% 
%     legend('Raw Data','Outliers','Least Squares Fit','Robust Fit')
%     xlabel('Area (mm^2)')
%     ylabel('Resistance (\Omega)')
%     grid on
% 
%     xlim([0 max(areas(:))])
%     ylim([floor(Z_min/100)*100 Z_max])
% end


%% Setup Video Figure

makeVideo = true;
filename_new_video = strcat(filename_video_raw, '_Matlab.mp4');

dpi_scale = 96/144; % matlab assumes 96dpi on Windows, but Dell M3800 is actually 144dpi
% figWidth = 1280*dpi_scale;
% figHeight = 720*dpi_scale;
figWidth = 1920*dpi_scale;
figHeight = 1080*dpi_scale;
gap = figHeight/100;

vid_raw.CurrentTime = start_time; % go to first frame

hFig = figure('units','pixels','Position',[100 10 figWidth figHeight],'MenuBar','none');

% video frame
hVid   = axes('units','pixels', 'Position',[0, 0, floor(vid_raw.Width*figHeight/vid_raw.Height), figHeight]);
hImage = imshow(vid_raw.readFrame,'border','tight', 'Parent',hVid);
hold on
hPoly = gobjects(4,1);
for ii=1:4
    hPoly(ii) = fill(hVid, polygons{1}{ii}(:,2), polygons{1}{ii}(:,1), colors(ii,:), 'FaceAlpha',0.5, 'EdgeColor',colors(ii,:));
%     hPoly(ii) = fill(hVid, polygons{1}{ii}(:,2), polygons{1}{ii}(:,1), cmap(Z_icolor(1,ii), :), 'FaceAlpha',0.5, 'EdgeColor',colors(ii,:));
end

% Z vs time
hZ = axes('units','pixels');
hold(hZ, 'on')
Z_times = (1:length(Z_sync))' / vid_raw.FrameRate; % [s]
ZaxisLimits = [0, Z_times(end), min(Z_sync(:))-50, max(Z_sync(:))+100];
axis(hZ, ZaxisLimits);
hZ_line   = gobjects(4,1);
hZ_marker = gobjects(4,1);
for ii=1:4
    hZ_line(ii) = animatedline(Z_times(1), Z_sync(1,ii), 'Color',colors(ii,:));
    hZ_marker(ii) = line(hZ, Z_times(1), Z_sync(1,ii), 'Color',colors(ii,:), 'LineStyle','none', 'Marker','o', 'MarkerSize',12);
end
xlabel(hZ, 'Time [s]');
ylabel(hZ, 'Z [\Omega]')
grid on

% Z vs Area
hArea_line   = gobjects(4,1);
hArea_marker = gobjects(4,1);
for ii=1:4
    hArea{ii}  = axes('units','pixels');
    hold(hArea{ii}, 'on')
    AreaAxisLimits = [0, 1.05*max(areas(:)), ZaxisLimits(3:4)];
    axis(hArea{ii}, AreaAxisLimits);
    hArea_line(ii) = animatedline(hArea{ii}, areas(1,ii), Z_sync(1,ii), 'Color', 'r', 'LineStyle','none', 'Marker','x', 'MarkerSize',5);
    hArea_marker(ii) = line(hArea{ii}, areas(1,ii), Z_sync(1,ii), 'Color', 'g', 'LineStyle','none', 'Marker','.', 'MarkerSize',30);
    ylabel(hArea{ii}, 'Z [\Omega]');
    xlabel(hArea{ii}, 'Area between electrodes [mm^2]');
    xticks(0:0.5:max(areas(:)))
    grid on
end

% position each plot
margins = max([hZ.TightInset; hArea{1}.TightInset]); % pick largest
hArea_width  = (figWidth-hVid.Position(3))/2 - margins(1)-margins(3) - gap;
hArea_height = (figHeight - 3*margins(2) - 3*margins(4))/3-gap;
hArea{1}.Position = [            hVid.Position(3)+  margins(1)+gap/2, hArea_height + 2*margins(2)+gap, hArea_width, hArea_height];
hArea{2}.Position = [hArea_width+hVid.Position(3)+2*margins(1)+gap,   hArea_height + 2*margins(2)+gap, hArea_width, hArea_height];
hArea{3}.Position = [            hVid.Position(3)+  margins(1)+gap/2,                  margins(2)+gap, hArea_width, hArea_height];
hArea{4}.Position = [hArea_width+hVid.Position(3)+2*margins(1)+gap,                    margins(2)+gap, hArea_width, hArea_height];
hZ.Position    = [hVid.Position(3)+margins(1)+gap/2, hArea{1}.Position(2)+hArea{1}.Position(4)+margins(2)+gap, figWidth-hVid.Position(3)-margins(1)-margins(3)-gap, hArea_height];

% set up VideoWriter
if makeVideo
    videoFWriter = vision.VideoFileWriter(fullfile(directory_name, filename_new_video), 'FileFormat','MPEG4', 'FrameRate',vid_raw.FrameRate, 'Quality',80);
end


%% plot/create video

vid_raw.CurrentTime = start_time; % go to first frame
currentFrame = 0; % reset count

tic;
while (hasFrame(vid_raw) && (currentFrame < size(areas,1)))

    % read next frame
    currentFrame = currentFrame + 1;
    vidFrame = readFrame(vid_raw);


    % Update Plots

    % current frame image
    hImage.CData = vidFrame; % faster than calling imshow

    for ii=1:4
        % Time vs Z plot
        addpoints(hZ_line(ii),   Z_times(currentFrame), Z_sync(currentFrame,ii));
        hZ_marker(ii).XData = Z_times(currentFrame);
        hZ_marker(ii).YData = Z_sync(currentFrame,ii);

        % polygon overlays
        hPoly(ii).XData = [polygons{currentFrame}{ii}(:,2); polygons{currentFrame}{ii}(1,2)]; % append first value to close polygon
        hPoly(ii).YData = [polygons{currentFrame}{ii}(:,1); polygons{currentFrame}{ii}(1,1)];
%         hPoly(ii).FaceColor = cmap(Z_icolor(currentFrame,ii), :);

        % Area vs Z plots
        addpoints(hArea_line(ii), areas(currentFrame,ii), Z_sync(currentFrame,ii));
        hArea_marker(ii).XData = areas(currentFrame,ii);
        hArea_marker(ii).YData = Z_sync(currentFrame,ii);
    end


    if makeVideo
        drawnow
        frame = getframe(hFig);
        videoFWriter(frame.cdata);
    else
        drawnow
%         pause(1/vid_raw.FrameRate);
    end

end

if makeVideo
    release(videoFWriter);
end
toc