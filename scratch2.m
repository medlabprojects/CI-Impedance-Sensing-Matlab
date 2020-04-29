%% 
% vid.CurrentTime = 15;
% im = vid.readFrame;
% cpuTime = timeit(@()segmentPixelsByColor(im, [1 0 0]), 1)
% gpuTime = timeit(@()segmentPixelsByColorGPU(im, [1 0 0]), 1)

%%
syms r a c1 c2 c3
% eqn1 = r == c1*a^c2 + c3;
eqn1 = r == c1*a^c2;
X = solve(eqn1, a, 'Real', true, 'IgnoreAnalyticConstraints',true)

%%

c(1) = exp(4.0986);
c(2) = -1.3791;
c(3) = 2638.7;

d(1) = exp(4.2);
d(2) = c(2);
d(3) = c(3);

A = linspace(0.2,1, 500);

Rc = c(1)*A.^c(2) + c(3);
Rd = d(1)*A.^d(2) + d(3);


figure(15); clf(15); hold on
plot(A, Rc, 'g', 'LineWidth', 2)
plot(A, Rd, 'r')
legend('Rc','Rd')
xlabel('Area (mm^2)')
ylabel('Resistance (\Omega)')
hold off

xlim([0 max(A)])
ylim([floor( min(c(3),d(3))/100 )*100 max([Rc Rd])])

%%

c(1) = exp(0);
c(2) = -1.5;
c(3) = 0.1;

d(1) = exp(0.5);
d(2) = c(2);
d(3) = c(3);

e(1) = exp(1);
e(2) = c(2);
e(3) = c(3);

A = linspace(0.05,1, 500);

Rc = c(1)*A.^c(2) + c(3);
Rd = d(1)*A.^d(2) + d(3);
Re = e(1)*A.^e(2) + e(3);


figure(16); clf(16); hold on
plot(A, Rc, 'g', 'LineWidth', 1.5)
plot(A, Rd, 'm', 'LineWidth', 1.5)
plot(A, Re, 'r', 'LineWidth', 1.5)

legend('Rc','Rd','Re')
xlabel('Area (mm^2)')
ylabel('Resistance (\Omega)')
hold off

xlim([0 max(A)])
ylim([floor( min(c(3),d(3))/100 )*100 max([Rc Rd])])



%% plot/create new video frames with data included

makeVideo = true;
filename_new_video = strcat(filename_video, '_Matlab.mp4');

Z = Z_sync;

% set up figure and axes
dpi_scale = 96/144; % matlab assumes 96dpi on Windows, but Dell M3800 is actually 144dpi
% figWidth = 1280*dpi_scale;
% figHeight = 720*dpi_scale;
figWidth = 1920*dpi_scale;
figHeight = 1080*dpi_scale;
gap = 5;
hFig = figure('units','pixels','Position',[100 10 figWidth figHeight],'MenuBar','none');

% video frame
hVid   = axes('units','pixels', 'Position',[0, 0, floor(vid.Width*figHeight/vid.Height), figHeight]);
hImage = imshow(vid.readFrame,'border','tight', 'Parent',hVid);

% Z vs time
hZ = axes('units','pixels');
hold(hZ, 'on')
Z_times = (1:length(Z)) / vid.FrameRate; % [s]
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
hArea.Position = [hVid.Position(3)+margins(1)+gap, margins(2)+gap,             figWidth-hVid.Position(3)-margins(1)-margins(3)-gap, figHeight/2-margins(2)-margins(4)-2*gap];
hZ.Position    = [hVid.Position(3)+margins(1)+gap, figHeight/2+margins(2)+gap, figWidth-hVid.Position(3)-margins(1)-margins(3)-gap, figHeight/2-margins(2)-margins(4)-2*gap];

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
