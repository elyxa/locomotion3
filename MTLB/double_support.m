function DS = double_support(rIC, lIC, rTC, lTC, KINtime)

if length(lTC) < length(rTC)
    num_steps = length(lTC);
else
    num_steps = length(rTC);
end
% (Right) gait cycle time
rIC_time = KINtime(rIC,2);
rGCT = [];
for i=1:num_steps
    rGCT(i) = rIC_time(i+1) - rIC_time(i);
end


for k=1:num_steps
% Terminal double stance

TDS(k) = (KINtime(rTC(k)) - KINtime(lIC(k)))./rGCT(k);

% Initial double stance
IDS(k) = (KINtime(lTC(k)) - KINtime(rIC(k)))./rGCT(k);

end
% Double stance
DS = TDS' + IDS';

end
