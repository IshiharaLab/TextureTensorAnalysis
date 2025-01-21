clear

path(path, 'code');

disp('=============')
disp('3.0 Load  data (.mat); e.g. "pos_link_step2.mat" ')

load pos_link_step2.mat

disp('=============')
disp('3.1 Example of Figures" ')
disp('I.  M and ROIs on first and last frames" ')
makeFig1(data)

disp('II.  M and ROIs on first and last frames" ')
plot_TrQandDT(data)

%%
disp('III.  The kinematic equation of texture tensor M1')
plot_dMelement(data)
