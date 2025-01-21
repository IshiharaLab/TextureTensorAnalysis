clear

path(path, 'code');

disp('=============')
disp('2.0 Load  data (.mat)')

[baseFileName, folder] = uigetfile('*.*', 'Specify an image sequence');
[pathstr,name,ext] = fileparts(baseFileName);
str0=[folder,name,ext];
data=load_step1data(str0);

%%%%%%%%%%%%%%%%%%% 
% Please define the area to be analyzed. 
% However, please be careful not to cut off the edges at the first and last frames.
% e.g.  h = impoly;
%        xy = getPosition(h);
%%%%%%%%%%%%%%%%%%% 
load input/xyROI.mat;
data.xyROI = xyROI;

%%%%%%%%%%%%%%%%%%% 
% Please fix any tracking errors if necessary.
%%%%%%%%%%%%%%%%%%% 

disp('=============')
disp('2.1 Division check')

data=division_cell(data);

%%
disp('=============')
disp('2.2 Separation of ROIs')

L =100;
data=sepa_ROI(data,L);

%%
disp('=============')
disp('2.3 Calculation of texture tensor')

data=cal_Texture(data);

%%
disp('=============')
disp('2.4 Save "analysis data"')

ffnam = strcat('pos_link_step2.mat');
[filename, PathName]  = uiputfile(ffnam,'Save file name (.mat)');
if filename == 0
    return
end
    
h=msgbox('Now saving ...','modal');
m = matfile(filename);
save(filename,'data')
close(h);
