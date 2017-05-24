function [stride_features, stridefeatures_names] = compute_kinematics(KIN, KINtime, IC, side)

kinematic_features = [];
kinfeatures_names = [];

n_steps = length(IC);

%% JOINTS - PLOTTING
%Relative movement between two contiguous segments
%Global frame through which local frames (limb) are expressed

% figure; 
% hold on;
% epochs=[];
% for i=1:n_steps-1
%     interval = [IC(i):IC(i+1)-1];
%     window = KIN.Ang.(side).(joint)(interval);
%     plot(window);
% end


%% Joint angles matrix
joint_ang_mat = [];
joint_ang_speed_mat =[];
joints ={ 'Hip';
          'Knee'; 
          'Ankle'; 
          'AbsAnkl'; 
          'Shoulder'; 
          'Elbow'; 
          'Wrist'; 
          'Neck'; 
          'Spine'; 
          'Head'; 
          'Thorax'; 
          'Pelvis'; 
          'FootProgress'};
for i=1:length(joints)
    joint_ang_mat =  [joint_ang_mat; KIN.Ang.(side).(joints{i})];
    joint_ang_speed_mat = [ joint_ang_speed_mat; diff(KIN.Ang.(side).(joints{i})) ];
end

%% velocity calcul
signal = KIN.Ang.(side).(joints{1});
pwelch(signal,[],[],[],100)
Nf = 50; 
Fpass = 10; 
Fstop = 12;

d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',100);

fvtool(d,'MagnitudeDisplay','zero-phase','Fs',100)

dt = KINtime(2)-KINtime(1);
vsignal = filter(d,signal)/dt;

tt = KINtime(1:end-delay);
vd = vsignal;
vd(1:delay) = [];

tt(1:delay) = [];
vd(1:delay) = [];

[pkp,lcp] = findpeaks(signal);
zcp = zeros(size(lcp));

[pkm,lcm] = findpeaks(-signal);
zcm = zeros(size(lcm));

subplot(2,1,1)
plot(KINtime,signal,KINtime([lcp lcm]),[pkp -pkm],'or')
xlabel('Time (s)')
ylabel('Displacement (cm)')
grid

subplot(2,1,2)
plot(tt,vd,KINtime([lcp lcm]),[zcp zcm],'or')
xlabel('Time (s)')
ylabel('Speed (cm/s)')
grid

adrift = filter(d,vdrift)/dt;

at = t(1:end-2*delay);
ad = adrift;
ad(1:2*delay) = [];

at(1:2*delay) = [];
ad(1:2*delay) = [];

subplot(2,1,1)
plot(tt,vd)
xlabel('Time (s)')
ylabel('Speed (cm/s)')
grid

subplot(2,1,2)
plot(at,ad)
ax = gca;
ax.YLim = 2000*[-1 1];
xlabel('Time (s)')
ylabel('Acceleration (cm/s^2)')
grid
%% Joints displacement matrix
joint_disp_mat = [];


% 
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

% planar covariation of elevation agles (ex: femur vs tibia)
%upper-lowe % left-right limb coordination
% ROM range of motion of knee and ankle joints

% joint moment and power, 



    % Intersegmental coordination
        % Coupling between adjacent segments
        % Linear coupling
    % Velocity of the joint angles
        % Amplitude of speed
        % Maximum speed
        % Minimum speed
    % Similarity of joint angles with pre-lesions
    
    % Joint angles
    
    % Limb oscillations
        % Maximum (backward position)
        % Minimum (backward position)
    
    % Trunk oscillations in the sagittal and frontal plane
    
    % Pelvis vertical oscillations
    
    % Limb endpoint trajectories
    
    % Timing of gait
    
    % Interlimb coordination
    
    % Variability
end 