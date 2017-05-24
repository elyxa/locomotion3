function [features, features_names] = com_energy(KIN, m) 


g = 9.80665;
%" the center of the pelvis, defined as the centroid of the triangle from the left anterior superior iliac spine, 
%the right anterior superior iliac spine, and the mid-point of the two posterior superior iliac spines, 
%could approximate the COM during walking ".
%Cite: 
%Eames MHA, Cosgrove A, Baker R. Comparing methods of estimating the total body centre of mass
%in three-dimensions in normal and pathological gaits. Human Movement Science. 1999;18:637?646

midPSI = (KIN.Pos.L.PSI + KIN.Pos.R.PSI)./2;
centroid = (midPSI + KIN.Pos.L.ASI + KIN.Pos.R.ASI)./3;

Vy = diff(centroid(1,:));
Vx = diff(centroid(2,:));
Vz = diff(centroid(3,:));

Ek_x = 0.5.*m.*Vx.^2
Ek_y = 0.5.*m.*Vy.^2
Ek_z = 0.5.*m.*(Vz.^2)

Ek = Ek_y + Ek_z + Ek_x;
Ep = m*g.*Vz;

Etot = Ek + Ep;

figure
plot(zscore(Ek));
hold on
plot(zscore(Ep));



features= [Ek; Ep];
features_names = ['Kinetic energy  '; 
                  'Potential energy'];
end