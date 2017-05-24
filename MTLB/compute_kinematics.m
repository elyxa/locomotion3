function [stride_features, stridefeatures_names] = compute_kinematics(KIN, KINtime, IC, side)


% THIS FUNCTION IS NOT REALLY USEFUL - JUST RANDOM STUFF :P


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