%% Import Data

% load trials 4/6/7, bin info, and trained network
load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\Dissertation data\Zsense_dissertation-data_trained-net.mat');
load('D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2020-01-23\Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6.mat', 'polygons');

% all are 4x1 cells with each cell 1xN
A = trial6.A;
Z = trial6.Z;
A_bin = trial6.A_bin;
A_bin_pred = trial6.A_bin_pred;


%% Import Videos

filename_video = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Impedance Sensing\Experiments\2020-01-23\Zsense_2020-01-27_phantom_Flex24_5x20_EA-1-5_trial6_crop.mp4';
vid     = VideoReader(filename_video);
start_time = 0.01;  % [s]
end_time   = vid.Duration; % [s]


%% Map predictions to colors

idx_pred_color = cell(4, 1); % each element specifies which row of class_colors to use
pred_color = cell(4,1); % CData arrays for Z vs A scatter plots

% find the corresponding color index for each measurement
for i_ch = 1:length(A_bin_pred)
    pred_color{i_ch} = zeros(length(A_bin{i_ch}), 3);
    for i_class = 1:length(bin_names)       
        % find all samples with this prediction
        idx = find(A_bin_pred{i_ch} == bin_names(i_class));
        idx_pred_color{i_ch}(idx) = i_class;
        pred_color{i_ch}(idx,:) = repmat(class_colors(i_class,:), [length(idx),1]);
    end
end


%% Setup Video Figure

makeVideo = false;
filename_new_video = strcat(filename_video, '_Matlab-Prediction.mp4');

dpi_scale = 96/144; % matlab assumes 96dpi on Windows, but Dell M3800 is actually 144dpi
figWidth = 1920*dpi_scale;
figHeight = 1080*dpi_scale;

vid.CurrentTime = start_time; % go to first frame

colors = distinguishable_colors(4);
hFig = figure('units','pixels','Position',[100 10 figWidth figHeight],'MenuBar','none');

% video frame
hVid   = axes('units','pixels', 'Position',[0, 0, floor(vid.Width*figHeight/vid.Height), figHeight]);
hImage = imshow(vid.readFrame,'border','tight', 'Parent',hVid);
hold on
hPoly = gobjects(4,1);
for ii=1:4
    % color by predicted location
    hPoly(ii) = fill(hVid, polygons{1}{ii}(:,2), polygons{1}{ii}(:,1), class_colors(idx_pred_color{ii}(1,ii),:), 'FaceAlpha',0.5, 'EdgeColor',colors(ii,:), 'LineWidth',1.2);
end

% Z vs time
hZ = axes('units','pixels');
hold(hZ, 'on')
Z_times = (1:length(Z{1}))' / vid.FrameRate; % [s]
ZaxisLimits = [0, Z_times(end), min(min(cell2mat(Z)))-50, max(max(cell2mat(Z)))+100];
axis(hZ, ZaxisLimits);
hZ_line   = gobjects(4,1);
hZ_marker = gobjects(4,1);
for ii=1:4
    hZ_line(ii) = animatedline(Z_times(1), Z{ii}(1), 'Color',colors(ii,:));
    hZ_marker(ii) = line(hZ,   Z_times(1), Z{ii}(1), 'Color',colors(ii,:), 'LineStyle','none', 'Marker','o', 'MarkerSize',12);
end
xlabel(hZ, 'Time [s]');
ylabel(hZ, 'R_a [\Omega]')
grid on

% Z vs Area
hArea_scatter = gobjects(4,1);
hArea_marker  = gobjects(4,1);
for ii=1:4
    hArea{ii}  = axes('units','pixels');
    hold(hArea{ii}, 'on')
    AreaAxisLimits = [0, 1.05*max(max(cell2mat(A))), ZaxisLimits(3:4)];
    axis(hArea{ii}, AreaAxisLimits);
    hArea_scatter(ii) = scatter(hArea{ii}, A{ii}, Z{ii}, '.', 'SizeData',5, 'CData', pred_color{ii});
    hArea_marker(ii) = line(hArea{ii}, A{ii}(1), Z{ii}(1), 'Color',colors(ii,:), 'LineStyle','none', 'Marker','.', 'MarkerSize',30);
    if (ii==1)||(ii==3)
        ylabel(hArea{ii}, 'R_a [\Omega]');
    else
        yticklabels([]);
    end
    xticks(0.5:0.5:AreaAxisLimits(2))
    if ii>2
        xlabel(hArea{ii}, 'Area [mm^2]');
    else
        xticklabels([]);
    end
    
    for i_class = 1:3
%         xline(hArea{ii}, bin_edges(i_class+1));
        fill(hArea{ii}, [bin_edges(i_class) bin_edges(i_class+1) bin_edges(i_class+1) bin_edges(i_class)],...
                        [AreaAxisLimits(3)  AreaAxisLimits(3)    AreaAxisLimits(4)    AreaAxisLimits(4)],...
             pred_colors(i_class,:), 'FaceAlpha',0.1, 'EdgeColor','none')
    end
    grid on
end

% position each plot
gap = 1.5 * figHeight/100; % percent of figure height
margins = max([hZ.TightInset; hArea{1}.TightInset]); % pick largest
hArea_width  = (figWidth-hVid.Position(3)- margins(1)- gap)/2;
hArea_height = (figHeight - 2*margins(2) - 2*margins(4))/3 - gap;
hArea{1}.Position = [            hVid.Position(3)+margins(1),     hArea_height + margins(2)+gap, hArea_width, hArea_height];
hArea{2}.Position = [hArea_width+hVid.Position(3)+margins(1)+gap/2, hArea_height + margins(2)+gap, hArea_width, hArea_height];
hArea{3}.Position = [            hVid.Position(3)+margins(1),                    margins(2)+gap/2, hArea_width, hArea_height];
hArea{4}.Position = [hArea_width+hVid.Position(3)+margins(1)+gap/2,                margins(2)+gap/2, hArea_width, hArea_height];
hZ.Position    = [hVid.Position(3)+margins(1), hArea{1}.Position(2)+hArea{1}.Position(4)+margins(2)+gap/2, figWidth-hVid.Position(3)-margins(1)-margins(3)-gap, hArea_height+gap];

% set up VideoWriter
if makeVideo
    videoFWriter = vision.VideoFileWriter(fullfile(directory_name, filename_new_video), 'FileFormat','MPEG4', 'FrameRate',vid.FrameRate, 'Quality',80);
end


%% plot/create video

update_rate = 1/30; % [s]
vid.CurrentTime = start_time; % go to first frame
currentFrame = 0; % reset count

a = tic;
while (hasFrame(vid) && (currentFrame < length(A{1})))

    % read next frame
    currentFrame = currentFrame + 1;
    vidFrame = readFrame(vid);


    % Update Plots

    % current frame image
    hImage.CData = vidFrame; % faster than calling imshow

    for ii=1:4
        % Time vs Z plot
        addpoints(hZ_line(ii),   Z_times(currentFrame), Z{ii}(currentFrame));
        hZ_marker(ii).XData = Z_times(currentFrame);
        hZ_marker(ii).YData = Z{ii}(currentFrame);

        % polygon overlays
        hPoly(ii).XData = [polygons{currentFrame}{ii}(:,2); polygons{currentFrame}{ii}(1,2)]; % append first value to close polygon
        hPoly(ii).YData = [polygons{currentFrame}{ii}(:,1); polygons{currentFrame}{ii}(1,1)];
        hPoly(ii).FaceColor = pred_colors(idx_pred_color(currentFrame,ii), :); % set color based on prediction

        % Area vs Z plots
        hArea_marker(ii).XData = A{ii}(currentFrame);
        hArea_marker(ii).YData = Z{ii}(currentFrame);
    end


    if makeVideo
        drawnow
        frame = getframe(hFig);
        videoFWriter(frame.cdata);
    else
        b = toc(a); % check timer
        if b > update_rate
            drawnow % update screen every X seconds
            a = tic; % reset timer after updating
        end
    end

end

if makeVideo
    release(videoFWriter);
end
toc