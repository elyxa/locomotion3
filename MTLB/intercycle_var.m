function PCI = intercycle_var(rIC, lIC)

num_steps = length(rIC);
phi=[];
for k=1:num_steps-1
   phi(k) = (rIC(k)-lIC(k))*360/(lIC(k+1)-lIC(k));
end

CV = std(phi) / mean(phi); 
PCI = mean(phi-180)*100/180+CV;

end