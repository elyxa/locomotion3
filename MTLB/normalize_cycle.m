function features = normalize_cycle(signal, IC)

n_steps = length(IC)
features=[];
figure;
hold on; 
for i=1:n_steps-1
    interval = [IC(i):IC(i+1)-1];
    if length(interval) > 100
        window = signal(interval);
        window = window(1:length(window)/100:end);
        plot(window);
        features = [features; window];
    end
 end

end