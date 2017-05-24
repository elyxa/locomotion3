function DS = double_support(rIC, lIC, rTC, lTC, KINtime)

% (Right) gait cycle time
rIC_time = KINtime(rIC,2);
rGCT = [];
for i=1:num_steps-1
    rGCT(i) = rIC_time(i+1) - rIC_time(i);
end

% Terminal double stance
TDS(k) = (rTC(k) - lIC(k))./rGCT(k);

% Initial double stance
IDS(k) = (lTC(k) - rIC(k))./rGCT(k);

% Double stance
DS = TDS + IDS;

end
