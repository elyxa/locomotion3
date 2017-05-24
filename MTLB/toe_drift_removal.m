function new_toe_signal = toe_drift_removal(KIN, KINtime, side)

%% drift removal for TOE - left
% this is very similar to "drift_removal"

t = KINtime(:,2);
y =  KIN.Pos.(side).TOE(3,:)';
[p, S] = polyfit(t,y,1);

y2 = polyval(p,t);
y3 = y-y2;
new_toe_signal = y3+(abs(y3(1)-y(1)));

% figure;
% plot(t,new_toe_signal);
% hold on
% plot(t,y)
end