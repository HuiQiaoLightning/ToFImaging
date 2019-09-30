% save depthmap to txt

addpath('..\functions')
sel_data = '0001';

fn = sprintf('results\\depth_%s',sel_data);
load(fn);
% depthmap, objdepth, grayvalue, tau_step, scenename);
% SPEEDOFLIGHT = 299792458e-6; % mm/ns
% k = SPEEDOFLIGHT*tau_step/2; % mm/depth
% d0 = 3600; % synchronizing signal delay offset
% distances = objdepth*k-d0;
distances = TOF2Distance(objdepth,tau_step);

%% show depthmap
% hdrgray = CompressHDR(grayvalue(1:100,1:163),0.01);
hdrgray = CompressHDR(grayvalue(:,1:163),0.008);
% hdrgray = grayvalue/max(grayvalue(:));
figure(1); image(hdrgray*255); 
% figure(1); imagesc(grayvalue); 
colormap(gray(256));colorbar

% dep1 = 0; dep2 = 1.8;
% objdepth = depthmap*tau_step;
% objdepth = min(objdepth,dep2);
% objdepth = max(objdepth,dep1);

% thr = max(grayvalue(:))*0.01;
% mask = grayvalue>thr;
% objdepth = objdepth.*mask + max(objdepth(:))*(~mask);
figure(2); imagesc(distances); 
colormap(jet(256));colorbar

fprintf('depthmap range %d ~ %d.\n',...
    min(depthmap(:)),max(depthmap(:)));

%% save files
fn = sprintf('results\\scene_%s_%s_depth.txt',sel_data,scenename); 
SaveArrayText(fn,distances,'%f');
% fn = sprintf('results\\scene_%s_%s_depth.png',sel_data,scenename); 
% imwrite(uint8(depthmap/max(depthmap(:))*256),fn,'bitdepth',8);
fn = sprintf('results\\scene_%s_%s_gray.png',sel_data,scenename); 
% imwrite(uint8(grayvalue/max(grayvalue(:))*256),fn,'bitdepth',8);
imwrite(uint8(hdrgray*255),fn,'bitdepth',8);
