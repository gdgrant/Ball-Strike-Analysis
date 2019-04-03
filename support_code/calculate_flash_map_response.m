function [flash_map_response,time_window, y, cp, cp_time_windows] = calculate_flash_map_response(session, get_single_flash_responses, calculate_choice_probabilities, green_or_cyan,time_window)

sorted_window = get_sorted_window(session.x);
response_window = [100 400]; 
if ~exist('time_window','var') || isempty(time_window)
    time_window = [-200 799];
end
trial_length = diff(time_window);
num_neurons = length(session.x(1).spike_times);
flash_map_response = zeros(num_neurons, 13, 26, trial_length, 'single');
flash_map_count = zeros(num_neurons, 13,26);
num_trials = length(session.x);
num_neurons = length(session.x(1).spike_times);

if get_single_flash_responses
    y = initialize_response_struct(num_neurons, trial_length);
else
    y = [];
end

for i = 1:num_trials
    if strcmp(session.x(i).trial_type, 'flash_map')
        [flash_time, flash_color, flash_zone, flash_pos, flash_color_name] = convert_flash_sequence(session.x(i), green_or_cyan);
        if length(flash_time) > 1 % need at least 2 flashes
            flash_pos(:,1) = flash_pos(:,1) + 1;
            flash_pos(:,2) = flash_pos(:,2) + 13.5;
            
            for j = 2:length(flash_time) % discard the first flash
                
                % make sure there's no fixation break release within 400 ms
                % of the flash and make sure we're not counting flashes
                % after the target flash
                if (session.x(i).TrialError ~= 3 || session.x(i).fixation_break - flash_time(j) > response_window(end)) ...
                        & (isempty(target_flash_time) || flash_time(j) < target_flash_time + 80)
                    
                    for k = 1:num_neurons
                        if i>=sorted_window(k,1) & i<=sorted_window(k,2)
                            st = session.x(i).spike_times{k} - flash_time(j) - time_window(1);
                            st = st(st>0 & st<=trial_length);
                            flash_map_response(k,flash_pos(j,1),flash_pos(j,2),st) = flash_map_response(k,flash_pos(j,1),flash_pos(j,2),st) + 1;
                            flash_map_count(k, flash_pos(j,1),flash_pos(j,2)) = flash_map_count(k,flash_pos(j,1),flash_pos(j,2)) + 1;
                            
                            if get_single_flash_responses
                                fn = [flash_color_name{j} num2str(flash_zone(j))];
                                y(k).(fn)(end+1,st) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

flash_map_count = reshape(flash_map_count,[num_neurons 13 26 1]);
flash_map_response = single(flash_map_response./repmat(flash_map_count,[1 1 1 trial_length]));

if calculate_choice_probabilities
    cp_time_windows = [[1:20:trial_length-100]' [101:20:trial_length]'];
    cp = calc_choice_probs(y, cp_time_windows);
else
    cp = [];
    cp_time_windows = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cp = calc_choice_probs(y, time_windows)

fn = fieldnames(y);
n = length(fn);
num_neurons = length(y);
num_time_windows = size(time_windows, 1);
cp = 0.5*ones(num_neurons, n, n, num_time_windows);

for t = 1:num_time_windows
    for g1 = 1:n
        for g2 = setdiff(1:n, g1)
            group1 = fn{g1};
            group2 = fn{g2};

            for k = 1:num_neurons
                resp1 =  mean(y(k).(group1)(:,time_windows(t,1):time_windows(t,2)),2);
                resp2 = mean(y(k).(group2)(:,time_windows(t,1):time_windows(t,2)),2);
                if size(resp1, 1) >=3 && size(resp2, 1) >=3
                    AUC = roc_masse(resp1, resp2);
                    cp(k, g1, g2, t) = AUC;
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initialize_response_struct(num_neurons, trial_length)

y = [];
for i = 1:num_neurons
    y(i).white1 = zeros(0,trial_length,'uint8');
    y(i).white2 = zeros(0,trial_length,'uint8');
    y(i).white3 = zeros(0,trial_length,'uint8');
end