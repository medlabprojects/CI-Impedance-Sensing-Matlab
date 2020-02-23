%%
tic
i_ch = 4;
updatedNet = net;
for ii=1:length(XValidation{i_ch})
    [updatedNet,label,score] = classifyAndUpdateState(updatedNet, XValidation{i_ch}(:,ii));
    labels(ii) = label;
end
toc

%%
% load('Zsensing_colormap.mat'); % load cmap (100 color levels)
% 
% % normalize each channel [0->1] based on global min/max
% Z_norm = (Z_sync - Z_min + eps) ./ (Z_max-Z_min);
% 
% % find the corresponding color index for each measurement
% Z_icolor = ceil(100*Z_norm);

%%
% % f1 = @() ZsensingSegmentFrame(RGB);
% % f2 = @() ZsensingSegmentFrameGPU(RGB);
% % t1 = timeit(f1, 3)
% % % t2 = timeit(f2, 3)
% % t3 = gputimeit(f2, 3)
% 
% %%
% % f1 = @() rgb2hsv(gpuArray(RGB));
% % f2 = @() gpuArray(rgb2hsv(RGB));
% % t1 = gputimeit(f1)
% % t2 = gputimeit(f2)
% 
% 
% %%
% vid.CurrentTime = 0.02;
% clear regions
% regions(length(frames)) = struct();
% for ii = 1:length(frames)
% 
%     if ~vid.hasFrame % first ensure we still have frames to read
%         break;
%     end
% 
%     % Segment current frame to find centers of electrode markers and outlines of modiolus/EA masks
%     RGB = vid.readFrame;
%     [regions(ii).polygons, regions(ii).areas] = ZsensingFindRegions(RGB);
% end
% 
% 
% 
% %%
% figure; hold on
% 
% scatter(A_all,Z_all) % Create scatter plot of resistance vs area
% % A_range = linspace(min(A_all),max(A_all),1000);
% A_range = linspace(0.1,2.5,1000);
% plot(A_range, F(x,A_range))
% % scatter(A_all(~idx_inliers), Z_all(~idx_inliers), 'xr') % mark outliers
% legend('data','fit')
% xlabel('Area (mm^2)')
% ylabel('Resistance (ohms)')
% 
% 
% %%
% modelfun = @(b,x) b(1).*x.^b(2) + b(3);
% beta0 = [250 -1.1 ];
% 
% %%
% ind = (1:length(Z_raw))';
% ind_sat = ind;
% Z_sat = Z_raw(:,4);
% Z_sat(Z_raw(:,4)<2000) = [];
% ind_sat(Z_raw(:,4)<2000) = [];
% Z_sat = interp1(ind_sat, Z_sat, ind);
% Z_smooth_sat = smooth(Z_sat, smooth_span, 'sgolay',3);
% figure; hold on;
% plot(Z_smooth_sat,'r')
% plot(Z_raw(:,4),'.k')
% plot(Z_sat,'.b')
% 
% 
% %%
% A = smooth(areas_raw(:,3), 30, 'sgolay', 4)';
% step = 8;
% A2x = 1:step:length(A);
% A2 = interp1(1:length(A), A, A2x);
% figure; hold on; grid minor;
% plot(diff(diff(A)), 'k')
% plot(A2x(2:end-1), diff(diff(A2)/step)/step,'b')

%%
% z=Z_trial4(:,3)-2000;
% z(z>5000) = 0;
% 
% a=areas(:,3);
% az = 300*a.^-1.2;
% 
% [C,lag] = xcorr(az-mean(az), z-mean(z), 1000);
% [~,I] = max(abs(C));
% z_shift = lag(I);
% 
% 
% figure
% ax(1) = subplot(1,2,1); 
% plot(lag,C,'k')
% ylabel('Amplitude')
% grid on
% title('Cross-correlation')
% 
% subplot(1,2,2);
% hold on
% plot(az,'g')
% if z_shift > 0
%     plot(z(z_shift:end), 'r')
% else
%     plot(z(1:end+z_shift), 'r')
% end

