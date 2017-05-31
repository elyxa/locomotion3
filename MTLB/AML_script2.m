clear all
close all
clc

% Select the c3d file to read:
cd('/Users/elisabettamessina/Desktop/Neuroscience_locomotion/MTLB/');
%addpath(['./']) %fold I'm in
addpath(['../Healthy/Healthy1']); %go up folder
addpath(['../Healthy/Healthy2']); %go up folder
addpath(['../SCI_HCU/SCI_HCU-noRobot2Crutches']);
addpath(['../SCI_HCU/SCI_HCU-RobotnoCrutches']); %go up folder
addpath(['../Healthy/Healthy2']);
addpath(['../Healthy/Student3and4']);


%% LOAD - SELECT
% Healthy subjects
healthy1 = cellstr(['Healthy1_Walk03.c3d';
                 'Healthy1_Walk04.c3d';
                 'Healthy1_Walk05.c3d';
                 'Healthy1_Walk06.c3d']);

healthy2 = cellstr(['Healthy2_Walk01.c3d';
                 'Healthy2_Walk03.c3d';
                 'Healthy2_Walk04.c3d';
                 'Healthy2_Walk05.c3d';
                 'Healthy2_Walk06.c3d' ]);
healthy3 = cellstr(['Healthy_Student3_normalspeed_01.c3d';
                 'Healthy_Student3_normalspeed_02.c3d';
                 'Healthy_Student3_normalspeed_03.c3d']);  
             
healthy4 = cellstr(['Healthy_Student4_normalspeed_01.c3d';
                    'Healthy_Student4_normalspeed_02.c3d';
                    'Healthy_Student4_normalspeed_03.c3d']);        
             
             
% Unhealthy subject - NO robot, 2 crutches
unhealthy_noAD = cellstr(['SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05.c3d';
                 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_06.c3d';
                 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_07.c3d';
                 'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_08.c3d' ]);     
gaitFile_noAD = cellstr(['SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05_GAIT.csv';
                  'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_06_GAIT.csv';
                  'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_07_GAIT.csv';
                  'SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_08_GAIT.csv']);

% Unhealthy subject - Robot, NO crutches
unhealthy_AD = cellstr(['SCI_HCU_20150505_02OVGa_AD_01.c3d';
                        'SCI_HCU_20150505_02OVGa_AD_02.c3d';    
                        'SCI_HCU_20150505_02OVGa_AD_03.c3d']);
gaitFile_AD = cellstr(['SCI_HCU_20150505_02OVGa_AD_01_GAIT.c3d';
                       'SCI_HCU_20150505_02OVGa_AD_02_GAIT.c3d';
                       'SCI_HCU_20150505_02OVGa_AD_03_GAIT.c3d']);



                   
                   
%% Loop for the healthy subject
f1_R=[];
f1_L=[];
f2=[];
f3 =[];
f4_L =[];
f4_R =[];
lAngles = [];
rAngles = [];
healthy_features = [];

for trial=1:length(healthy1) %CHANGE WITH PATIENT
    
    Filename = char(healthy1(trial)); %CHANGE WITH PATIENT
    [KIN,KINfreq,EMG,EMGfreq,Measures, FirstFrame]=load_c3d(Filename);
    
    % TIME (column #1 = sample# , column #2 = time (s))
    KINtime = [1:length(KIN.Pos.L.ANK(1,:))]' +FirstFrame;
    KINtime = [KINtime, (KINtime-1)./KINfreq];
    EMGtime = [1 :length(EMG.L.TA(1,:))]' + (FirstFrame*EMGfreq/KINfreq);
    EMGtime = [EMGtime, (EMGtime-1)./EMGfreq];

  
    % Find Heel-strikes
    [~, lHS] = findpeaks(-KIN.Pos.L.HEE(3,:), 'MinPeakDistance',50);
    [~, rHS] = findpeaks(-KIN.Pos.R.HEE(3,:), 'MinPeakDistance',50);
    
    % Find Toe-off
    lTO = find_toeoff(KIN.Pos.L.TOE(3,:), lHS);
    rTO = find_toeoff(KIN.Pos.R.TOE(3,:), rHS);
    
    % Fix Toe Off
    st_r=find(rHS(rHS<min(rTO)));
    rHS=rHS(st_r(end):end);
    st_l=find(lHS(lHS<min(lTO)));
    lHS=lHS(st_l(end):end);
      
    % Heel and toe signals without drift
    R_heel_zsignal = heel_drift_removal(KIN, rHS, 'R');
    L_heel_zsignal = heel_drift_removal(KIN, lHS, 'L');

    R_toe_zsignal = toe_drift_removal(KIN, KINtime, 'R');
    L_toe_zsignal = toe_drift_removal(KIN, KINtime, 'L');
   

%COMPUTE FEATURES    
   %Spatiotemporal features
    %f1_R =[f1_R; compute_spatiotemporal(KIN, KINtime, rHS, rTO, 'R')];
    %f1_L =[f1_L; compute_spatiotemporal(KIN, KINtime, lHS, lTO, 'L')];
    %f1_names = {'GCT' 'Cadence' 'StrideLength' 'GaitSpeed' 'Swing%'}; % 1,2,3,4,5
    %f2 = [f2; double_support(rHS, lHS, rTO, lTO, KINtime)]; %6
    %f2_name = {'Double support'};

    %Clearance
    %f4_R = [f4_R; compute_clearance(KIN, rHS, rTO,'R')]; 
    %f4_L = [f4_L; compute_clearance(KIN,lHS, lTO, 'L')]; 
    %f4_names = {'maxHeelClearance', 'maxToeClearance', 'minToeClearance'}; % 9, 10, 11
    
    %Elevation angles
     %lAngles =[lAngles; elevation_angles(KIN, lHS, 'L')];
     %rAngles = [rAngles; elevation_angles(KIN, rHS, 'R')];
     % Now we need to choose the features - remove rAngles to avoid
     % recomputing the features for the same cycles multiple times
  
    f3 = [f3; com_energy(KIN, Measures.Mass, rHS)]; %there is something wrong with kE
    f3_names = {'Kinetic energy','Potential energy'}; %7, 8
   
end

%%
count=size(f1_L,1)+size(f1_R,1);
Cycle = [];
Subject = [];
for i=1:count
Cycle = [Cycle; {['Healthy gait cycle' num2str(i)]}];
Subject = [Subject; {['Healthy']}];
end

[writeR{1:size(f1_R,1)}] = deal('R');
[writeL{1:size(f1_L,1)}] = deal('L');

Side = [writeR, writeL]';

f1 = [f1_R; f1_L];

T = table(Cycle, Subject, Side, Condition, f1(:,1), f1(:,2), f1(:,3), f1(:,4), f1(:,5)); % Trial_name),  f3, f4, f5);
filen = 'patientdata.xlsx';
writetable(T,filen,'Sheet',1,'Range','A1');



%% Loop for the UNhealthy subject
f1_R=[];
f1_L=[];
f2=[];
f3 =[];
f4_L =[];
f4_R =[];
lAngles = [];
rAngles = []; 
unhealthy_features = [];
for trial=1:length(unhealthy_noAD) %CHANGE WITH PATIENT
    
    Filename = char(unhealthy_noAD(trial)); %CHANGE WITH PATIENT
    [KIN,KINfreq,EMG,EMGfreq,Measures, FirstFrame]=load_c3d(Filename);
    gaitFile = char(gaitFile_noAD(trial)); %CHANGE WITH PATIENT
    gaitFile = readtext(gaitFile);
    
    % TIME (column #1 = sample# , column #2 = time (s))
    KINtime = [1:length(KIN.Pos.L.ANK(1,:))]' +FirstFrame;
    KINtime = [KINtime, (KINtime-1)./KINfreq];
    EMGtime = [1 :length(EMG.L.TA(1,:))]' + (FirstFrame*EMGfreq/KINfreq);
    EMGtime = [EMGtime, (EMGtime-1)./EMGfreq];

  
    % Find gait events
   [rHS, lHS, rTO, lTO ] = read_gait_events(gaitFile, KINtime, KINfreq, FirstFrame);

    % Heel and toe signals without drift
    R_heel_zsignal = heel_drift_removal(KIN, rHS, 'R');
    L_heel_zsignal = heel_drift_removal(KIN, lHS, 'L');

    R_toe_zsignal = toe_drift_removal(KIN, KINtime, 'R');
    L_toe_zsignal = toe_drift_removal(KIN, KINtime, 'L');
    
    %COMPUTE FEATURES    
   %Spatiotemporal features
    f1_R =[f1_R; compute_spatiotemporal(KIN, KINtime, rHS, rTO, 'R')];
    f1_L =[f1_L; compute_spatiotemporal(KIN, KINtime, lHS, lTO, 'L')];
    f1_names = {'GCT' 'Cadence' 'StrideLength' 'GaitSpeed' 'Swing%'}; % 1,2,3,4,5
    f2 = [f2; double_support(rHS, lHS, rTO, lTO, KINtime)]; %6
    f2_name = {'Double support'};
   %Energies
    f3 = [f3; com_energy(KIN, Measures.Mass, rHS)]; %there is something wrong with kE
    f3_names = {'Kinetic energy','Potential energy'}; %7, 8
    
   %Clearance
    %f4_R = [f4_R; compute_clearance(KIN, rHS, rTO,'R')]; 
    %f4_L = [f4_L; compute_clearance(KIN,lHS, lTO, 'L')]; 
    f4_names = {'maxHeelClearance', 'maxToeClearance', 'minToeClearance'}; % 9, 10, 11
    
    %Elevation angles
%     lAngles =[lAngles; elevation_angles(KIN, lHS, 'L')];
%     rAngles = [rAngles; elevation_angles(KIN, rHS, 'R')];
     % Now we need to choose the features
    
    unhealthy_features = [f1_R, f1_L, f2, f3, f4_R, f4_L];
end













% Plotting:
%% Plot left heel marker (vertical vs longitudinal)
figure, plot(KIN.Pos.L.HEE(2,:),L_heel_zsignal,'linewidth',2);
title('Heel marker (vertical vs longitudinal)');
xlabel('Longitudinal pos');
ylabel('Vertical pos');

%% Plot left/right heel (+ toe) markers vertical vs time (+ heel-strikes)
figure, plot(KINtime(:,2), L_heel_zsignal,'linewidth',2);
hold on, plot(KINtime(:,2), L_toe_zsignal,'linewidth',2);
plot( [KINtime(lHS(1),2) KINtime(lHS(1),2)], [0 30], 'k--');
plot( [KINtime(lTO(1),2) KINtime(lTO(1),2)], [0 30],'m--');
title('Left heel and toe markers vertical vs time + HS and TO');
legend('Heel marker z pos', 'Toe marker z pos', 'Heel strike', 'Toe off');
%plot( [KINtime(lHS,2) KINtime(lHS,2)], [0 30], 'k--');
%plot( [KINtime(lTO,2) KINtime(lTO,2)], [0 30],'m--');
xlabel('Time');
ylabel('Vertical pos');


figure, plot(KINtime(:,2),  R_heel_zsignal,'linewidth',2);
hold on, plot(KINtime(:,2), R_toe_zsignal,'linewidth',2);
plot( [KINtime(rHS(1),2) KINtime(rHS(1),2)], [0 30], 'k--');
plot( [KINtime(rTO(1),2) KINtime(rTO(1),2)], [0 30],'m--');
title('Right heel and toe markers vertical vs time + HS and TO');
legend('Heel marker z pos', 'Toe marker z pos', 'Heel strike', 'Toe off');
plot( [KINtime(rHS,2) KINtime(rHS,2)], [0 30], 'k--');
plot( [KINtime(rTO,2) KINtime(rTO,2)], [0 30],'m--');
xlabel('Time');
ylabel('Vertical pos');




