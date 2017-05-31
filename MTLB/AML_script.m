
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

%% LOAD HEALTHY

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

%% LOAD UNHAELTHY + read the gait file to set gait events

% Unhealthy subject - Robot, NO crutches
    %Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05_GAIT.csv');

    %Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_06.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_06_GAIT.csv');

    %Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_07.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_07_GAIT.csv');
    
    %Filename = 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_08.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_08_GAIT.csv');


% Unhealthy subject - NO robot, 2 crutches
    %Filename = 'SCI_HCU_20150505_02OVGa_AD_01.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_02OVGa_AD_01_GAIT.c3d');

    %Filename = 'SCI_HCU_20150505_02OVGa_AD_02.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_02OVGa_AD_02_GAIT.c3d');

    %Filename = 'SCI_HCU_20150505_02OVGa_AD_03.c3d';
    %gaitFile = readtext('SCI_HCU_20150505_02OVGa_AD_03_GAIT.c3d');


%% Load the file
[KIN,KINfreq,EMG,EMGfreq,Measures, FirstFrame]=load_c3d(Filename);
%FirstFrame
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
KINtime = [1:length(KIN.Pos.L.ANK(1,:))]' +FirstFrame;
KINtime = [KINtime, (KINtime-1)./KINfreq];
EMGtime = [1 :length(EMG.L.TA(1,:))]' + (FirstFrame*EMGfreq/KINfreq);
EMGtime = [EMGtime, (EMGtime-1)./EMGfreq];
 

%% Plot left heel marker (vertical vs longitudinal)
figure, plot(KIN.Pos.L.HEE(2,:),KIN.Pos.L.HEE(3,:),'linewidth',2);

%% HEALTHY - detect gait events and show a graph

% Find Heel-strikes
[~, L_HEEL_STRIKE] = findpeaks(-KIN.Pos.L.HEE(3,:), 'MinPeakDistance',50);
[~, R_HEEL_STRIKE] = findpeaks(-KIN.Pos.R.HEE(3,:), 'MinPeakDistance',50);

% Find Toe-off

[~, R_Toe_STRIKE] = findpeaks(KIN.Pos.R.TOE(3,:),'MinPeakDistance',70);
 
for i = 1:(length(R_Toe_STRIKE)-1)
    n=1;
    p=1;
    for j = round((R_Toe_STRIKE(i+1)- R_Toe_STRIKE(i))/2+R_Toe_STRIKE(i)):R_Toe_STRIKE(i+1)
    if KIN.Pos.R.TOE(3,j+1)>0.4*KIN.Pos.R.TOE(3,R_Toe_STRIKE(i+1))
         if KIN.Pos.R.TOE(3,j)<0.4*KIN.Pos.R.TOE(3,R_Toe_STRIKE(i+1))
            Toe_off(p) = n;  
            p=p+1;
         end
    end
     n=n+1;
    end 
    R_TOE_OFF(i) = Toe_off(1)+round((R_Toe_STRIKE(i+1)- R_Toe_STRIKE(i))/2+R_Toe_STRIKE(i));
end
% [~, Toe_off] = findpeaks(-KIN.Pos.R.TOE(3,R_Toe_STRIKE(i):R_Toe_STRIKE(i+1)),'MinPeakDistance',R_Toe_STRIKE(i+1)-R_Toe_STRIKE(i)-1);
% R_TOE_OFF(i) = Toe_off+R_Toe_STRIKE(i);
% end
% Take many points and take the max
[~, L_Toe_STRIKE] = findpeaks(KIN.Pos.L.TOE(3,:),'MinPeakDistance',70);
 
for i = 1:(length(L_Toe_STRIKE)-1)
    n=1;
    p=1;
    Toe_off = 0;
    for j = round((L_Toe_STRIKE(i+1)- L_Toe_STRIKE(i))/2+L_Toe_STRIKE(i)):L_Toe_STRIKE(i+1)
    if KIN.Pos.L.TOE(3,j+1)>0.4*KIN.Pos.L.TOE(3,L_Toe_STRIKE(i+1))
         if KIN.Pos.L.TOE(3,j)<0.4*KIN.Pos.L.TOE(3,L_Toe_STRIKE(i+1))
            Toe_off(p) = n;  
            p=p+1;
         end
    end
     n=n+1;
    end 
    L_TOE_OFF(i) = Toe_off(1)+round((L_Toe_STRIKE(i+1)- L_Toe_STRIKE(i))/2+L_Toe_STRIKE(i));
end

%FOR TOE OFF

R_HEEL_STRIKE=R_HEEL_STRIKE(R_HEEL_STRIKE<max(R_TOE_OFF));
R_TOE_OFF=R_TOE_OFF(R_TOE_OFF>min(R_HEEL_STRIKE));
 
st_r=find(R_HEEL_STRIKE(R_HEEL_STRIKE<min(R_TOE_OFF)));
R_HEEL_STRIKE=R_HEEL_STRIKE(st_r(end):end);
 
R_HEEL_STRIKE=R_HEEL_STRIKE(R_HEEL_STRIKE<max(R_TOE_OFF));
 
L_HEEL_STRIKE=L_HEEL_STRIKE(L_HEEL_STRIKE<max(L_TOE_OFF));
L_TOE_OFF=L_TOE_OFF(L_TOE_OFF>min(L_HEEL_STRIKE));
 
st_l=find(L_HEEL_STRIKE(L_HEEL_STRIKE<min(L_TOE_OFF)));
L_HEEL_STRIKE=L_HEEL_STRIKE(st_l(end):end);
 
L_HEEL_STRIKE=L_HEEL_STRIKE(L_HEEL_STRIKE<max(L_TOE_OFF));


%% UNHEALTHY - gait events

[R_HEEL_STRIKE, L_HEEL_STRIKE, R_TOE_OFF, L_TOE_OFF ] = read_gait_events(gaitFile, KINtime) ;

R_HEEL_STRIKE= R_HEEL_STRIKE - FirstFrame;
L_HEEL_STRIKE = L_HEEL_STRIKE - FirstFrame;
R_TOE_OFF = R_TOE_OFF - FirstFrame;
L_TOE_OFF = L_TOE_OFF - FirstFrame;
%% Undrift heel and toe
R_heel_signal = heel_drift_removal(KIN, R_HEEL_STRIKE, 'R');
L_heel_signal = heel_drift_removal(KIN, L_HEEL_STRIKE, 'L');

R_toe_signal = toe_drift_removal(KIN, KINtime, 'R');
L_toe_signal = toe_drift_removal(KIN, KINtime, 'L');


%% Plot left/right heel (+ toe) markers vertical vs time (+ heel-strikes)
figure, plot(KINtime(:,2), L_heel_signal,'linewidth',2);
hold on, plot(KINtime(:,2), L_toe_signal,'linewidth',2);
plot( KINtime(L_HEEL_STRIKE,2), 1 ,'o');
plot( KINtime(L_TOE_OFF,2), 1 ,'*');
title('Left heel and toe markers vertical vs time + HS and TO');

figure, plot(KINtime(:,2),  R_heel_signal,'linewidth',2);
hold on, plot(KINtime(:,2), R_toe_signal,'linewidth',2);
plot( KINtime(R_HEEL_STRIKE,2), 1 ,'o');
plot( KINtime(R_TOE_OFF,2), 1 ,'*');
title('Right heel and toe markers vertical vs time + HS and TO');


% The Heel Strike have been evaluated as the minima of the z component of
% the heel marker. It could also be visualized through the z component of the same marker 
% (when the plot becomes flat)
% The Toe-offs have been evaluated as...


%% compute gait parameters of your choice (at least 30 params based on kinematics and 5 based on
%EMGs), and explain them - load patient data
 

% compute the avarage temporal parameters - 
% change 'R' with 'L' to get left side

features_matrix = []; %problem: how to deal with the difference in length??

% We get 4 features x 3 (1 per gait cycle) :
[temp_features, temp_feat_names] = compute_spatiotemporal(KIN, KINtime, R_HEEL_STRIKE, 'R');
[temp_features, temp_feat_names] = compute_spatiotemporal(KIN, KINtime, R_HEEL_STRIKE, 'L');

% We get 3 features x 4 (1 per peak)
[clearence_features, clea_feat_names] = compute_clearance(KIN, R_HEEL_STRIKE, R_TOE_OFF,'R'); 
[clearence_features, clea_feat_names] = compute_clearance(KIN, L_HEEL_STRIKE, L_TOE_OFF, 'L'); 

% We get 1 feature
PCI = intercycle_var( R_HEEL_STRIKE,  L_HEEL_STRIKE); %phase coordination index (inter-cycle variability)

% We get 3 features x 414 (n of frames)
[elev_angles, elev_angles_names] = elevation_angles(KIN, 'R');
[elev_angles, elev_angles_names] = elevation_angles(KIN, 'L');

% We get 2 feature x 413 (n of frames -1)
[en_features, en_names] = com_energy(KIN, Measures.Mass);

% We get
DS = double_support(R_HEEL_STRIKE,  L_HEEL_STRIKE, R_TOE_OFF, L_TOE_OFF, KINtime);

%% Other stuff we still need
% each line - each gait cycle
 %normalize 0 100

% pelvis  obliquity
%trunk inclination and rigidity (frontal and horizontal planes)
% step width
% transverse foot angles 
% symmetry
% velocities
%Foot Angular velocity
% stance and swing vector





%% EMG: high pass filter (10-2k Hz)
% HP 30 Hz to remove artifact movements
% 50 Hz filter (you can do at the same time)
% Rectification (neg to pos)
% Low pass filter (denaturate the signal) - threshold to detect beginning
% and end of the signal
% beginning and end of the signal (duration)
% mean, integral, max


%% compute gait parameters
%(text?) file with patient and all the features to compute


%% save all parameters in a single Excel sheet. The final excel sheet that you should produce is
%constructed as follows: first columns = text columns with the conditions 
%(less than 10 columns in total); then columns of parameters; first line = header line, 
%other lines = 1 line per gait cycle

%% load the .txt file in GSTAT (you can find the executable on the moodle) - perform the PCA

%% discuss, extract relevant parameters, discuss, conclude


