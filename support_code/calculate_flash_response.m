function [flash_response, y, cp, cp_time_windows, flash_count, fn, LFP_power] = calculate_flash_response(session, desired_block_number, get_single_flash_responses, ...
	calculate_choice_probabilities, green_or_cyan, time_window, power_freq_bands, power_LFP_channels, trial_type)

	% Handle and set defaults as necessary
	if ~exist('get_single_flash_responses','var') || isempty(get_single_flash_responses)
		get_single_flash_responses = true;
	end
	if ~exist('calculate_choice_probabilities','var') || isempty(calculate_choice_probabilities)
		calculate_choice_probabilities = false;
	end

	if ~exist('power_freq_bands','var')
		power_freq_bands = [];
	end
	if ~exist('power_LFP_channels','var')
		power_LFP_channels = [];
	end

	if ~exist('trial_type','var') || isempty(trial_type)
		trial_type = 'ball_strikes';
	end
	if ~exist('green_or_cyan','var') || isempty(green_or_cyan)
		green_or_cyan = 'green';
	end
	if ~exist('desired_block_number','var') || isempty(desired_block_number)
		desired_block_number = 1;
	end

	if ~exist('time_window','var') || isempty(time_window)
		time_window = [-200 400];
	end

	if ~exist('power_freq_bands','var')
		power_freq_bands = [];
		power_LFP_channels = [];
		step_size = 4;
		num_freq_bands = length(power_freq_bands);
		num_ch = length(power_LFP_channels);
	else
		step_size = 1;
		num_freq_bands = [];
		num_ch = [];
	end

	% Get the sorted window
	sorted_window = get_sorted_window(session.x);

	% Set parameters
	response_window = [100 400];
	num_trials = length(session.x);
	num_neurons = length(session.x(1).spike_times);
	trial_length = diff(time_window);

	% Set up data storage
	flash_response = zeros(num_neurons, 3, 13, 26, trial_length, 'single');
	flash_count = zeros(num_neurons, 3,13,26);
	flash_count_neuron = zeros(num_neurons, 1);

	% Identify monkey-dependent parameters
	if strcmp(session.Monkey, 'Wahwah')
		num_dead_flashes = 1;
	elseif strcmp(session.Monkey, 'Quincy')
		num_dead_flashes = 1;
	end

	% Make the appropriate data deposition struct
	if get_single_flash_responses
		[y, ~, fn] = initialize_response_struct(num_neurons, trial_length, num_freq_bands, num_ch, step_size);
	else
		y = [];
		fn = [];
		LFP_step_size = 1;
	end

	% Start a flash counter, then iterate over the available trials
	flash_number = 0;
	for i = 1:num_trials

		% Isolate the trial
		trial = session.x(i);

		% If this trial is not appropriate to analyze, continue to the next one
		if ~strcmp(trial.trial_type, trial_type) || trial.block_number ~= desired_block_number
			continue
		end

		% Retrieve the flash sequence information
		% [flash_time, flash_color, flash_zone, flash_pos, ~, target_flash_time] = convert_flash_sequence(trial, green_or_cyan);
		flash_time = trial.flash_time;
		flash_color = trial.flash_color;
		flash_zone = trial.flash_zone;
		flash_pos = trial.flash_pos;
        
		target_flash_inds = find(flash_color == flash_zone & flash_zone < 3);
        if ~isempty(target_flash_inds)
            target_flash_time = flash_time(target_flash_inds(1));
        else
            target_flash_time = [];
        end
		
		% If selected, compute LFP power
		if ~isempty(power_freq_bands)
			LFP = zeros(length(power_LFP_channels), length(trial.LFP{1}));
			for j = 1:length(power_LFP_channels)
				LFP(j,:) = trial.LFP{power_LFP_channels(j)};
			end
			ws = complex_morlet_transform(LFP, power_freq_bands, 1000, 6);
		end

		% If there are an insufficient number of flashes, continue to the next trial
		if length(flash_time) < 2
			continue
		end
		
		% Translate the flash positions to index units
		flash_pos(:,1) = flash_pos(:,1) + 1;
		flash_pos(:,2) = flash_pos(:,2) + 13.5;

		% Iterate over relevant flashes
		for j = num_dead_flashes+1:length(flash_time)

			% Ensure that there's no fixation break before the flash, and that
			% flashes after the target flash are not used.  If either of these
			% conditions fail, continue to the next flash
            % NOTE: response_window(1)   ---->   response_window(2)
            % NOTE: < target_flash_time + 80   ---->   <= target_flash_time
            if ~((trial.TrialError ~= 3 || trial.fixation_break - flash_time(j) > response_window(1)) ...
               && (isempty(target_flash_time) || flash_time(j) < target_flash_time + 80) ...
               && (isempty(trial.lever_up) || trial.lever_up - flash_time(j) > response_window(1)))
                continue
            end
			
			% Iterate the flash number (only occurs if the flash is valid)
			flash_number = flash_number + 1;

			% Determine the lever_up state for this flash
			if ~isempty(trial.lever_up) && (trial.lever_up - flash_time(j) > response_window(1) ...
					& trial.lever_up - flash_time(j) <= response_window(end))
				lever_up = 1;
			else
				lever_up = 0;
			end
			
			% Determine the fixation_break state for this flash
			if ~isempty(trial.fixation_break) && (trial.fixation_break - flash_time(j) > response_window(1) ...
					& trial.fixation_break - flash_time(j) <= response_window(end))
				fix_break = 1;
			else
				fix_break = 0;
			end
			
			% Determine the target_flash state for this flash
			if flash_color(j)<=2 & flash_color(j) == flash_zone(j)
				target_flash = 1;
			else
				target_flash = 0;
			end
			
			% Determine the LFP_power for this flash
			if ~isempty(power_freq_bands) & get_single_flash_responses
				u = [flash_time(j)+time_window(1):flash_time(j)+time_window(end)];
				u = u(LFP_step_size:LFP_step_size:end); 
				LFP_power(flash_number, :, :, :) = squeeze(ws(:,u,1,:));
			end
			
			% Iterate over the available neurons to record this information
			for n = 1:num_neurons

				% Only record this trial if it is within the sorted window
				if i>=sorted_window(n,1) && i<=sorted_window(n,2)

					% Iterate the neuron's flash count
					flash_count_neuron(n) = flash_count_neuron(n) + 1;
					flash_count(n, flash_color(j), flash_pos(j,1), flash_pos(j,2)) = ...
						flash_count(n, flash_color(j), flash_pos(j,1), flash_pos(j,2)) + 1;

					% Calculate the neuron's spike times and flash response
					st = trial.spike_times{n} - flash_time(j) - time_window(1);
					st = st(st>0 & st<=trial_length);
					flash_response(n, flash_color(j), flash_pos(j,1), flash_pos(j,2), st) = ...
						flash_response(n, flash_color(j), flash_pos(j,1), flash_pos(j,2), st) + 1;

					% If we want the single flash responses, record into the output struct
					if get_single_flash_responses
						y(n).spike_resp(flash_count_neuron(n), st) = true;
						y(n).flash_pos(flash_count_neuron(n), :)   = flash_pos(j,:);
						y(n).flash_col(flash_count_neuron(n))      = flash_color(j);
						y(n).flash_zone(flash_count_neuron(n))     = flash_zone(j);
						y(n).lever_up(flash_count_neuron(n))       = lever_up;
						y(n).fix_break(flash_count_neuron(n))      = fix_break;
						y(n).trial_error(flash_count_neuron(n))    = trial.TrialError;
						y(n).target_flash(flash_count_neuron(n))   = target_flash;
						y(n).flash_number(flash_count_neuron(n))   = flash_number;

						if strcmp(trial_type, 'flash_map')
							y(n).trial_type(flash_count_neuron(n)) = 1;
						elseif strcmp(trial_type, 'ball_strikes')
							y(n).trial_type(flash_count_neuron(n)) = 2;
						end

					end
				end
			end
		end
	end

	% Reformulate the flash count and response
	flash_count = reshape(flash_count,[num_neurons 3 13 26 1]) + 1e-6;
	flash_response = single(flash_response./repmat(flash_count,[1 1 1 1 trial_length]));

	% Further plug output data into the output struct for each neuron
	if get_single_flash_responses

		% Iterate over neurons
		for n = 1:num_neurons
			y(n).spike_resp    = y(n).spike_resp(1:flash_count_neuron(n),:);
			y(n).flash_pos     = y(n).flash_pos(1:flash_count_neuron(n),:);
			y(n).flash_col     = y(n).flash_col(1:flash_count_neuron(n));
			y(n).flash_zone    = y(n).flash_zone(1:flash_count_neuron(n));
			y(n).lever_up      = y(n).lever_up(1:flash_count_neuron(n));
			y(n).fix_break     = y(n).fix_break(1:flash_count_neuron(n));
			y(n).flash_number  = y(n).flash_number(1:flash_count_neuron(n));
			y(n).target_flash  = y(n).target_flash(1:flash_count_neuron(n));
			y(n).trial_error   = y(n).trial_error(1:flash_count_neuron(n));
			y(n).trial_type    = y(n).trial_type(1:flash_count_neuron(n));

			% If appropriate, format LFP data
			if step_size > 1
				temp = conv2(single(y(n).spike_resp), ones(1,LFP_step_size),'same');
				y(n).spike_resp = uint8(temp(:,LFP_step_size:LFP_step_size:end));
			end
		end

		% If power_freq_bands has been filled, operate upon it accordingly
		if ~isempty(power_freq_bands)
		   LFP_power = LFP_power(1:flash_number,:,:,:);
		   % Make the zscore and average across the channels
		   for j = 1:length(power_LFP_channels)
			  LFP_temp  = reshape(LFP_power(:,j,:,:),[size(LFP_power,1)*size(LFP_power,3) length(power_freq_bands)]);
			  LFP_temp = zscore(LFP_temp);
			  LFP_power(:,j,:,:) = reshape(LFP_temp,[size(LFP_power,1) size(LFP_power,3) length(power_freq_bands)]);
		   end
		   LFP_power = angle(squeeze(mean(LFP_power, 2)));
		end

	end

	% Calculate the choice probabilities if desired, otherwise return empties
	if calculate_choice_probabilities
		cp_time_windows = [[1:20:trial_length-100]' [101:20:trial_length]'];
		cp = calc_choice_probs(y, cp_time_windows);
	else
		cp_time_windows = [];
		cp = [];
	end

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
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, LFP_power, fn] = initialize_response_struct(num_neurons, trial_length, num_freq_bands, num_channels, step_size)

	y = [];
	max_num_flashes = 15000;
	for i = 1:num_neurons
		y(i).spike_resp = zeros(max_num_flashes,trial_length,'uint8');
		y(i).spike_resp = logical(y(i).spike_resp);
		y(i).trial_error = zeros(max_num_flashes,1,'uint8');
		y(i).flash_pos = zeros(max_num_flashes,2,'uint8');
		y(i).flash_zone = zeros(max_num_flashes,1,'uint8');
		y(i).flash_col = zeros(max_num_flashes,1,'uint8');
		y(i).lever_up = zeros(max_num_flashes,1,'uint8');
		y(i).fix_break = zeros(max_num_flashes,1,'uint8');
		y(i).target_flash = zeros(max_num_flashes,1,'uint8');
		y(i).flash_number = zeros(max_num_flashes,1,'single');
		y(i).trial_type = zeros(max_num_flashes,1,'uint8');
		if ~isempty(num_freq_bands)
			LFP_power = zeros(max_num_flashes,num_channels, floor(trial_length/step_size),num_freq_bands,'single');
		else
			LFP_power = [];
		end
	end

	fn = fieldnames(y);
end