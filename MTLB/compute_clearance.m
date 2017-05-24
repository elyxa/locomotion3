function [features, features_names] = compute_clearance(KIN, side)

% FOOT CLEARANCE

% Max heel clearance
%figure; plot(KIN.Pos.(side).HEE(3,:))
[max_HC, max_HC_ind] = findpeaks(KIN.Pos.(side).HEE(3,:), 'MinPeakDistance',70);
%hold on; plot(max_HC_ind, max_HC, 'o');

% Max toe (undrifted)
undrifted_signal = drift_removal(KIN.Pos.(side).TOE(3,:));

%figure; plot(undrifted_signal);
%hold on; plot((KIN.Pos.(side).TOE(3,:)));
[max_TC, max_TC_ind] = findpeaks(undrifted_signal, 'MinPeakDistance',70);
%hold on; plot(max_TC_ind, max_TC, 'o');

% Min toe (undrifted)
a = smooth(undrifted_signal, 10);
%figure; plot(a);
[min_TC, min_TC_ind] = findpeaks(-a,  'MaxPeakWidth',20, 'MinPeakWidth',10);
%hold on; plot(min_TC_ind, -min_TC, 'o');
min_TC = -min_TC';

features= [max_HC; max_TC; min_TC];
features_names = ['maxHeelClearance'; 'maxToeClearance '; 'minToeClearance '];

end