
function new_heel_signal = heel_drift_removal(KIN, IC, side)

% remove drift for HEE 
%Needs IC points
y1=KIN.Pos.(side).HEE(3,:);
query = KIN.Pos.(side).HEE(3,IC);
line = linspace(query(1),query(end),414);
y2=KIN.Pos.(side).HEE(3,:)-line;
%diff=(y1(1)-y2(1));
new_heel_signal = y2+abs(y2(1)-y1(1));

% figure;
% plot(new_heel_signal);
% hold on
% plot(y1)

end