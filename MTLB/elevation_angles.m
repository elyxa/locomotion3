function [elev_angles, names] = elevation_angles(KIN, side)

% ELEVATION ANGLES

%Elevation angles are defined as the oscillation of the limb with respect to 
%the direction of gravity. In comparison to join angles, elevation angle are 
%remarkably similar across a range of speeds

y_hip = KIN.Pos.(side).ASI(3,:);
x_hip = KIN.Pos.(side).ASI(2,:);
y_knee = KIN.Pos.(side).KNE(3,:);
x_knee = KIN.Pos.(side).KNE(2,:);
y_ankle = KIN.Pos.(side).ANK(3,:);
x_ankle = KIN.Pos.(side).ANK(2,:);
y_toe = KIN.Pos.(side).TOE(3,:);
x_toe = KIN.Pos.(side).TOE(2,:);

%pelvis ? ????

%thigh
a = x_knee - x_hip;
b = y_hip - y_knee;
theta_thigh = atan(a./b);

%shank
 
a = x_knee - x_ankle ;
b = y_knee - y_ankle;
theta_shank = - atan(a./b);

%foot - not working
a = x_toe - x_ankle;
b = y_ankle - y_toe;
theta_foot = atan(a./b);

a = find(theta_foot<0)
theta_foot(a) = theta_foot(a) + pi;

figure; hold on;
subplot(3,1,1), plot(theta_foot);
title('Foot elevation angle');
subplot(3,1,2), plot(theta_shank);
title('Shank elevation angle');
subplot(3,1,3), plot(theta_thigh);
title('Thigh elevation angle');
    elev_angles = [theta_thigh; theta_shank];
    names = ['Thigh elev. angle'; 'Shank elev. angle'];
end