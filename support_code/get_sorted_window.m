function sorted_window = get_sorted_window(x)

num_neurons = length(x(1).spike_times);
sorted_window = zeros(num_neurons, 2);
num_trials = length(x);
for i = 1:num_neurons
    for j = 1:num_trials
        if ~isempty(x(j).spike_times{i})
            sorted_window(i,1) = j;
            break
        end
    end
    for j = num_trials:-1:1
        if ~isempty(x(j).spike_times{i})
            sorted_window(i,2) = j;
            break
        end
    end
end