path(path, 'code');

disp('=============')
disp('1.0 Load image')
    
% % Browse for the image file.
% [baseFileName, folder] = uigetfile('*.*', 'select the skeletonized image data (.tiff)');
% [pathstr,name,ext] = fileparts(baseFileName);
% fullImageFileName = fullfile(folder, baseFileName);
%fullImageFileName = '/Users/namba/Desktop/upload_code/input/Skeleton-Cregion-WT3.tif';
fullImageFileName = 'input/Skeleton-Cregion-WT3.tif';

% disp('=============')
% disp('1.5 Load PIV data')
% 
% [baseFileName, folder] = uigetfile('*.*', 'select the "PIVlab.mat"');
% PIVname=[folder,baseFileName];

%PIVname='/Users/namba/Desktop/upload_code/input/PIVlab.mat';
PIVname='input/PIVlab.mat';
data.flgPIV="on";
   
data = load_Image(data,fullImageFileName,PIVname);
