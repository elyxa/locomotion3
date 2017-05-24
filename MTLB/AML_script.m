clear all
close all
clc

% Select the c3d file to read:
cd('/Users/elisabettamessina/Desktop/Neuroscience_locomotion/MTLB/')
%addpath(['./']) %fold I'm in
addpath(['../Healthy/Healthy1']) %go up folder
addpath(['../Healthy/Healthy2']) %go up folder
addpath(['../SCI_HCU/SCI_HCU-noRobot2Crutches'])
addpath(['../SCI_HCU/SCI_HCU-RobotnoCrutches']) %go up folder
addpath(['../Healthy/Healthy2'])

%% LOAD - select only one of these at time =)

% Healthy subject #1
Filename ='Healthy1_Walk03.c3d';
%Filename ='Healthy1_Walk04.c3d';
%Filename ='Healthy1_Walk05.c3d';
%Filename ='Healthy1_Walk06.c3d';

% Healthy subject #2
%Filename ='Healthy2_Walk01.c3d';
%Filename ='Healthy2_Walk03.c3d';
%Filename ='Healthy2_Walk04.c3d';
%Filename ='Healthy2_Walk05.c3d';
%Filename ='Healthy2_Walk06.c3d';

% Unhealthy subject - Robot, NO crutches
%Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05.c3d';
%Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_06.c3d';
%Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_07.c3d';
%Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_08.c3d';


% Unhealthy subject - NO robot, 2 crutches
%Filename = 'SCI_HCU_20150505_02OVGa_AD_01.c3d';
%Filename = 'SCI_HCU_20150505_02OVGa_AD_02.c3d';
%Filename = 'SCI_HCU_20150505_02OVGa_AD_02.c3d';

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

%% read GAIT file to set the gait events 


% create for healthy - plot initial contact and last contact 


%% Plot left heel marker (vertical vs longitudinal)

figure, plot(KIN.Pos.L.HEE(2,:),KIN.Pos.L.HEE(3,:),'linewidth',2)


%% detect gait events (by eye) and show a graph


% Heel strikes
[~, L_HEEL_STRIKE] = findpeaks(-KIN.Pos.L.HEE(3,:), 'MinPeakDistance',50)
[~, R_HEEL_STRIKE] = findpeaks(-KIN.Pos.R.HEE(3,:), 'MinPeakDistance',50)

% Undrift heel and toe
R_heel_signal = heel_drift_removal(KIN, R_HEEL_STRIKE, 'R');
L_heel_signal = heel_drift_removal(KIN, L_HEEL_STRIKE, 'L');

R_toe_signal = toe_drift_removal(KIN, KINtime, 'R');
L_toe_signal = toe_drift_removal(KIN, KINtime, 'L');


%Note: SHOW TOE OFF GRAPHICALLY

% Plot left/right heel markers vertical vs time (+ heel-strikes)
figure, plot(KINtime(:,2), L_heel_signal,'linewidth',2);
hold on, plot(KINtime(:,2), L_toe_signal,'linewidth',2);
plot( KINtime(L_HEEL_STRIKE,2), 1 ,'o');

figure, plot(KINtime(:,2),  R_heel_signal,'linewidth',2);
hold on, plot(KINtime(:,2), R_toe_signal,'linewidth',2);
plot( KINtime(R_HEEL_STRIKE,2), 1 ,'o');


% The Heel Strike have been evaluated as the minima of the z component of
% the heel marker. It could also be visualized through the z component of the same marker 
% (when the plot becomes flat)


%% compute gait parameters of your choice (at least 30 params based on kinematics and 5 based on
%EMGs), and explain them - load patient data
 

% compute the avarage temporal parameters - 
% change 'R' with 'L' to get left side

features_matrix = []; %problem: how to deal with the difference in length??

% We get 4 features x 3 (1 per gait cycle) :
[temp_features, temp_feat_names] = compute_spatiotemporal(KIN, KINtime, R_HEEL_STRIKE, 'R');
[temp_features, temp_feat_names] = compute_spatiotemporal(KIN, KINtime, R_HEEL_STRIKE, 'L');

% We get 3 features x 4 (1 per peak)
[clearence_features, clea_feat_names] = compute_clearance(KIN, 'R'); 
[clearence_features, clea_feat_names] = compute_clearance(KIN, 'L'); 

% We get 1 feature
PCI = intercycle_var( R_HEEL_STRIKE,  L_HEEL_STRIKE); %phase coordination index (inter-cycle variability)

% We get 3 features x 414 (n of frames)
[elev_angles, elev_angles_names] = elevation_angles(KIN, 'R');
[elev_angles, elev_angles_names] = elevation_angles(KIN, 'L');

% We geet 2 feature x 413 (n of frames -1)
[en_features, en_names] = com_energy(KIN, Measures.Mass);



%% Other stuff we still need


% pelvis  obliquity
%trunk inclination and rigidity (frontal and horizontal planes)
% step width
% transverse foot angles 
% step period
% double stance - toe-off
% symmetry
% velocities
%Foot Angular velocity
% stance and swing vector
% joint moment and power




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


