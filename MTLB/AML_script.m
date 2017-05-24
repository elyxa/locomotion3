clear all
close all
clc

%% load healthy raw data

% Select the c3d file to read:
cd('/Users/elisabettamessina/Desktop/EPFL2sem/locomotion/Assignments/Neuroscience_locomotion/MTLB/')
%addpath(['./']) %fold I'm in
addpath(['../Healthy/Healthy1']) %go up folder
addpath(['../SCI_HCU/SCI_HCU-noRobot2Crutches'])
addpath(['../SCI_HCU/SCI_HCU-RobotnoCrutches']) %go up folder
addpath(['../Healthy/Healthy2'])
Filename = 'Healthy1_Walk03.c3d';

% Load the file
[KIN,KINfreq,EMG,EMGfreq,Measures]=load_c3d(Filename);

% Data you have (KIN=kinematic / EMG=Electromyographic):
% - KIN.Pos.X.YYY = 3D position of marker YYY on side X
%                   - X = L (left), R (right)
%                   - YYY = TOE (toe), HEE (heel), ANK (ankle), TIB (tibia,
%                   on shank), KNE (knee), THI (thigh), ASI (anterior
%                   superior iliac crest) or PSI (posterior superior iliac
%                   crest), ...
% - KIN.Ang.X.Yyyy = flexion/extension angle at joint center Yyyy on side X
%                   - X = L (left), R (right)
%                   - Yyyy = Knee or Ankle
% - KIN.JCs.X.Yyyy = 3D position of joint center Yyyy on side X
%                   - X = L (left), R (right)
%                   - Yyyy = HJC (Hip), KJC (Knee), AJC (Ankle) or FJC (Foot)
% - KINfreq = sampling frequency of KIN data
% - EMG.X.Yyy = EMG activity of muscle Yyy on side X
%                   - X = L (left), R (right)
%                   - Yyy = MG (Medial Gastrocnemius), TA (Tibialis Anterior)
% - EMGfreq = sampling frequency of EMG data
% - Measures = anthropometric data
%

% Axes are such that: X = lateral (sideward), Y = longitudinal (forward), Z = vertical (upward)


% TIME (column #1 = sample# , column #2 = time (s))
KINtime = [1:length(KIN.Pos.L.ANK(1,:))]';
KINtime = [KINtime, (KINtime-1)./KINfreq];
EMGtime = [1:length(EMG.L.TA(1,:))]';
EMGtime = [EMGtime, (EMGtime-1)./EMGfreq];
  
%% Plot left heel marker (vertical vs longitudinal)
figure, plot(KIN.Pos.L.HEE(2,:),KIN.Pos.L.HEE(3,:),'linewidth',2)

%plot3(KIN.Pos.L.HEE(1,:), KIN.Pos.L.HEE(2,:), KIN.Pos.L.HEE(3,:))


%% detect gait events (by eye) and show a graph

[~, L_HEEL_STRIKE] = findpeaks(-KIN.Pos.L.HEE(3,:), 'MinPeakDistance',50)
[~, R_HEEL_STRIKE] = findpeaks(-KIN.Pos.R.HEE(3,:), 'MinPeakDistance',50)

%Note: SHOW TOE OFF GRAPHICALLY

% Plot left/right heel markers vertical vs time (+ heel-strikes)
figure, plot(KINtime(:,2), KIN.Pos.L.HEE(3,:),'linewidth',2);
hold on, plot(KINtime(:,2), KIN.Pos.L.TOE(3,:),'linewidth',2);
plot( KINtime(L_HEEL_STRIKE,2), 1 ,'o');

figure, plot(KINtime(:,2), KIN.Pos.R.HEE(3,:),'linewidth',2);
hold on, plot(KINtime(:,2),KIN.Pos.R.TOE(3,:),'linewidth',2);
plot( KINtime(R_HEEL_STRIKE,2), 1 ,'o');


% The Heel Strike have been evaluated as the minima of the z component of
% the heel marker. It could also be visualized through the z component of the same marker 
% (when the plot becomes flat)

%% drift removal for TOE
% automatic
t = KINtime(:,2);
y =  KIN.Pos.L.TOE(3,:)';
[p, S] = polyfit(t,y,1);
plot(t,y)
title('Plot of y Versus t')

y2 = polyval(p,t);
figure
plot(t,y,'o',t,y2)
title('Plot of Data (Points) and Model (Line)')

y3 = y-y2;
figure;
plot(t,y3+(abs(y3(1)-y(1))))
hold on
plot(t,y)

%% random stuff
% figure
% plot(KIN.Pos.L.HEE(3,:));
% [t, y] = ginput(4);
% [p, S] = polyfit(t,y,1);
% y2 = polyval(p,t);
% plot(t,y,'o',t,y2)
%% remove drift for HEE
%Needs IC points
y1=KIN.Pos.L.HEE(3,:)
query = KIN.Pos.L.HEE(3,L_HEEL_STRIKE);
line = linspace(query(1),query(end),414);
y2=KIN.Pos.L.HEE(3,:)-line;
%diff=(y1(1)-y2(1));
figure;
plot(y2+abs(y2(1)-y1(1)));
hold on
plot(y1)
%% compute gait parameters of your choice (at least 30 params based on kinematics and 5 based on
%EMGs), and explain them - load patient data
 

% compute the avarage temporal parameters - 
% change 'R' with 'L' to get left side
[temp_features, temp_feat_names] = compute_spatiotemporal(KIN, KINtime, R_HEEL_STRIKE, 'R');

[clearence_features, clea_feat_names] = compute_clearance(KIN, 'R'); 

PCI = intercycle_var( R_HEEL_STRIKE,  L_HEEL_STRIKE); %phase coordination index (inter-cycle variability)

[elev_angles, elev_angles_names] = elevation_angles(KIN, 'R');

compute_kinematics(KIN, R_HEEL_STRIKE, 'R');
%[en_features, en_names] = com_energy(Measures, Mass);


%symmetry index (need swing times)
%turning angles

%potential energy - com_energy
%kinetic energy - com_energy
%%



%% read GAIT file to set the gait events (you can copy-paste in MTLB)

% ciao = readtable('SCI_HCU_20150505_02OVGa_AD_01_GAIT.csv');


% event = toe off
% create for healthy - plot initial contact and last contact 
% joint center - estimation
% use marker position and angle elevation

%% EMG: high pass filter (10-2k Hz)
% HP 30 Hz to remove artifact movements
% 50 Hz filter (you can do at the same time)
% Rectification (neg to pos)
% Low pass filter (denaturate the signal) - threshold to detect beginning
% and end of the signal
% beginning and end of the signal (duration)
% mean, integral, max
%% compute gait parameters
%file with patient and all the features to compute
%% save all parameters in a single Excel sheet. The final excel sheet that you should produce is
%constructed as follows: first columns = text columns with the conditions (less than 10 columns in total); then columns of parameters; first line = header line, other lines = 1 line per gait cycle
%% load the .txt file in GSTAT (you can find the executable on the moodle) - perform the PCA
%% discuss, extract relevant parameters, discuss, conclude


%% 