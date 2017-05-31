function TO = find_toeoff(toe_signal, IC)
    
TO=[];
n_steps = length(IC)-1
features=[];
for i=1:n_steps
    interval = [IC(i):IC(i+1)-1];
    if length(interval) > 100
        window = toe_signal(interval);
        a = diff(window);
        [~,col] = find(a>0.2)
        TO = [TO, col(1)+IC(i)-1]; %not sure about -1
    end
end
    

end
