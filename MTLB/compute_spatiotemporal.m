function [temporal_features, tempfeatures_names] = compute_spatiotemporal(KIN, KINtime, IC, TO, side)


temporal_features = [];
tempfeatures_names = [];

num_steps = length(IC);
IC_time = KINtime(IC,2);



% gait cycle
GCT=[];

for i=1:num_steps-1
    GCT(i) = IC_time(i+1) - IC_time(i);
end

% cadence (steps per min)
Cad = 120./GCT;

% % step period - need to provide both legs
% lStP = [];
% rStP = [];
% 
% for i=1:num_steps
%     lStP(i) = rIC_time(i) - lIC_time(i);
%     rStP(i) = rIC_time(i) - lIC_time(i);
% end


SL=[];
SV=[];
SW=[];
for i=1:num_steps-1
    
    %stride length
    SL(i) = KIN.Pos.(side).HEE(2, IC(i+1)) - KIN.Pos.(side).HEE(2, IC(i)) ;

    % gait speed (stride velocity) 
    SV(i) = SL(i)/GCT(i);
    
    % swing / (swing+stance)
    SW(i) = TO(i);
end




temporal_features = [GCT', Cad', SL', SV', SW'];
tempfeatures_names = {'GCT' 'Cadence' 'StrideLength' 'GaitSpeed' 'Swing%'};




end