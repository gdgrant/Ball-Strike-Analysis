function y = ball_strike_single_unit_stats_single_session_downsample(session, green_or_cyan)

	% Add the current path
	addpath('./')

	% Set some parameters
	num_reps = 10;
	percentiles = [80,85,90,95];

	% Set some defaults
	desired_block_number = 1;
	get_single_flash_responses = 1;
	calculate_choice_probabilities = 0;
	skip_flash_map_data = true;

	if ~exist('green_or_cyan','var') || isempty(green_or_cyan)
		green_or_cyan = 'green';
	end

	% Set up environment parameters and output struct
	num_neurons = length(session.x(1).spike_times);

	y = [];
	y.percentiles = percentiles;
	y.num_neurons = num_neurons;
	% y.time_window = [-200 400];
	% y.time_windows1 = [[-199:10:301]' [-100:10:400]'];
	y.time_window = [-100 300];
	y.time_windows1 = [[-99:10:201]' [0:10:300]'];
	% y.time_windows2 = [[-199:10:351]' [-150:10:400]'];

	options = define_options();
	roc_descriptions = define_roc_description();

	% Parse basic neuron properties
	y.session_date   = session.ExperimentDate;
	y.channel_num    = cat(1,session.NeuronInfo(:).Channel);
	y.iso_qual       = cat(1,session.NeuronInfo(:).IsolationQuality);
	y.depth          = cat(1,session.NeuronInfo(:).Depth);
	[y.waveform_width, y.waveform_amp, ~, ~] = analyze_waveforms(session);

	%%%%% Start analysis %%%%%

	% Extract spike and flash information for each task
	if skip_flash_map_data
		[y.ball_strikes_response, y.single_flash_resp, ~, ~, y.ball_strikes_flash_count] = calculate_flash_response(session, desired_block_number, ...
			get_single_flash_responses, calculate_choice_probabilities, green_or_cyan, y.time_window,[],[],'ball_strikes');
		y.single_trial_combined = y.single_flash_resp;
		disp('Passing through <y.single_flash_resp>, ignoring combine_single_trial_data().')
	else
		[y.ball_strikes_response, y.single_flash_resp, ~, ~, y.ball_strikes_flash_count] = calculate_flash_response(session, desired_block_number, ...
			get_single_flash_responses, calculate_choice_probabilities, green_or_cyan, y.time_window,[],[],'ball_strikes');
		[y.flash_map_response, y.single_flash_resp_flash_map,~,~,y.flash_map_flash_count] = calculate_flash_response(session, desired_block_number, ...
			get_single_flash_responses, calculate_choice_probabilities, green_or_cyan, y.time_window,[],[],'flash_map');
		y.single_trial_combined = combine_single_trial_data(y.single_flash_resp_flash_map, y.single_flash_resp);
	end

	% Calculate delayed memory saccade stats
	if false
		[y.mem_sac_response, y.mem_sac_pref_angle, y.mem_sac_p_val, y.mem_sac_mod_depth, y.mem_sac_epochs] = get_mem_saccade_response(session);
	end
		
	% Calculate correct/failed ROC values
	if false
		[y.green_dp2, y.red_dp2] = calculate_detect_prob(session, desired_block_number, y.time_windows2, green_or_cyan);
		[y.green_dp1, y.red_dp1] = calculate_detect_prob(session, desired_block_number, y.time_windows1, green_or_cyan);
	end

	% Get the response tally for each neuron
	for n=1:num_neurons
		y.num_flashes(n,:) = size(y.single_trial_combined(n).spike_resp, 1);
	end
		
	%%% Ball Strikes 1
	colors = 1:3;
	bs1 = calculate_single_neuron_stats(y.single_trial_combined, y.time_windows1, y.time_window(1), colors, num_reps, percentiles, options, roc_descriptions, green_or_cyan);
	% y.comparisons = bs1.comparisons;
	
	%%% Ball Strikes 2
	% bs2 = calculate_single_neuron_stats(y.single_flash_resp, y.time_windows2, y.time_window(1), colors, num_reps, options, roc_descriptions, green_or_cyan);
	
	%%% Flash Map 1
	% colors = 3;
	% fm1 = calculate_single_neuron_stats(y.single_trial_combined, y.time_windows1, y.time_window(1), colors, num_reps, percentiles, options, roc_descriptions, green_or_cyan);
	
	%%% White Comparison 1
	% colors = 1:2;
	% wc1 = calculate_single_neuron_stats(y.single_trial_combined, y.time_windows1, y.time_window(1), colors, num_reps, options, roc_descriptions, green_or_cyan);
	
	%%% Flash Map 2
	% fm2 = calculate_single_neuron_stats(y.single_flash_resp_flash_map, y.time_windows2(5:5:end,:), y.time_window(1), colors, num_reps, options, roc_descriptions, green_or_cyan);
	
	% Add all analysis data to the output struct by fieldname
	fn = fieldnames(bs1);
	for i = 1:length(fn)
		if ~strcmp(fn{i},'comparisons')
			new_fn1 = ['bs_' fn{i} '1'];
			y.(new_fn1) = bs1.(fn{i});
			% new_fn2 = ['bs_' fn{i} '2'];
			% y.(new_fn2) = bs2.(fn{i});
			% new_fn1 = ['fm_' fn{i} '1'];
			% y.(new_fn1) = fm1.(fn{i});
			% new_fn2 = ['fm_' fn{i} '2'];
			% y.(new_fn2) = fm2.(fn{i});
			% new_fn1 = ['white_' fn{i}];
			% y.(new_fn1) = wc1.(fn{i});
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function single_trial_combined = combine_single_trial_data(single_flash_resp_flash_map, single_flash_resp)

	single_trial_combined = [];
	num_neurons = length(single_flash_resp);
	for n = 1:num_neurons
		if length(single_flash_resp_flash_map(n).flash_zone) > 100 && length(single_flash_resp(n).flash_zone) > 100
			fn = fieldnames(single_flash_resp(n));
			for j = 1:length(fn)
				single_trial_combined(n).(fn{j}) = cat(1, single_flash_resp_flash_map(n).(fn{j}), single_flash_resp(n).(fn{j}));
			end
			n1 = size(single_flash_resp_flash_map(n).flash_col, 1);
			n2 = size(single_flash_resp(n).flash_col, 1);
			single_trial_combined(n).trial_type = cat(1,ones(n1,1), 2*ones(n2,1));
		else
			fn = fieldnames(single_flash_resp(n));
			for j = 1:length(fn)
				single_trial_combined(n).(fn{j}) = [];
			end
			single_trial_combined(n).trial_type = [];
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = calculate_single_neuron_stats(single_flash_resp, time_windows, align_time, colors, num_reps, percentiles, options, roc_descriptions, green_or_cyan)

	% Align time windows
	time_windows = time_windows - align_time;

	% Set up environment
	num_neurons = length(single_flash_resp);
	num_windows = size(time_windows, 1);
	num_pcts = length(percentiles);

	% Put together some parameters
	window_size = [3 3];
	num_comparisons = 357;
	% [x_pos, y_pos] = meshgrid(0:12,-12.5:12.5);

	% Establish output data struct
	y = [];
	y.RF_mask            = zeros(num_neurons, num_reps, 3, num_windows, 2, 13, 26, 'single');
	y.RF_ROC             = zeros(num_neurons, num_reps, 3, num_windows, 2, 'single');
	y.RF_MI              = zeros(num_neurons, num_reps, 3, num_windows, 2, 'single');
	y.RF_fit_params      = zeros(num_neurons, num_reps, 3, num_windows, 7, 'single');
	y.RF_gaussian_fit    = zeros(num_neurons, num_reps, 3, num_windows, 13, 26, 'single');
	y.RF_R2              = zeros(num_neurons, num_reps, 3, num_windows, 'single');
	y.flash_roc          = zeros(num_neurons, num_reps, num_comparisons, num_windows, 'single');
	y.perc_mask          = zeros(num_neurons, num_reps, 3, num_windows, num_pcts, 13, 26, 'logical');
	% y.smoothed_perc_mask = zeros(num_neurons, num_reps, 3, num_windows, num_pcts, 13, 26, 'logical');
	y.mean_spiking       = zeros(num_neurons, num_reps, 3, num_windows, 13, 26, 'single');
	y.mean_spiking_ds    = zeros(num_neurons, num_reps, 3, num_windows, 13, 26, 'single');
	y.flash_count        = zeros(num_neurons, num_reps, 3, num_windows, 13, 26, 'single');
	y.flash_count_ds     = zeros(num_neurons, num_reps, 3, num_windows, 13, 26, 'single');

	% Iterate over neurons and repetitions for analysis
	for n = 1:num_neurons
		for r = 1:num_reps
			
			% Downsample flashes relevant to the neuron
			disp(['Calculating stats for neuron ' num2str(n) ', repetition ' num2str(r) '.'])
			[spike_resp_ds, flash_col_ds, flash_pos_ds, flash_zone_ds, target_flash_ds] = ...
				downsample_target_flashes(single_flash_resp(n).spike_resp, single_flash_resp(n).flash_pos, ...
				single_flash_resp(n).flash_col, single_flash_resp(n).flash_zone, green_or_cyan);
			
			% If there are no sampled flashes, go to the next repetition or neuron
			if isempty(spike_resp_ds)
				continue
			end
			
			% Iterate over colors
			for col = colors

				% Isolate spikes with only the correct associated flash color
				ind = find(single_flash_resp(n).flash_col == col);
				spike_resp_current = single_flash_resp(n).spike_resp(ind,:);
				flash_pos_current  = single_flash_resp(n).flash_pos(ind,:);

				% Isolate downsampled spikes with only the correct associated flash color
				ind = find(flash_col_ds == col);
				spike_resp_ds_current = spike_resp_ds(ind,:);	% This was casted as a double at one point
				flash_pos_ds_current  = flash_pos_ds(ind,:);
				
				% Smooth data for gaussian fit
				% [mean_spike_count_ds, flash_count_ds_current] = ...
				% 	convert_single_trial_to_mean_resp(spike_resp_ds_current, flash_pos_ds_current);
				% [smoothed_spike_resp_ds_current, smoothed_flash_count_ds_current] = ...
				% 	smooth_data(mean_spike_count_ds, x_pos, y_pos, flash_count_ds_current);
				
				% Iterate over time windows
				for t = 1:num_windows
					
					% Sum up spikes in the current window
					timeframe = time_windows(t,1):time_windows(t,2);
					spike_count_current_ds = squeeze(sum(spike_resp_ds_current(:,timeframe),2))';
					spike_count_current = squeeze(sum(spike_resp_current(:,timeframe),2))';

					% smoothed_spike_count_current_ds = squeeze(sum(smoothed_spike_resp_ds_current(:,:,time_windows(t,1):time_windows(t,end)),3))';
					
					%%% Calculate mean spiking
					[y.mean_spiking(n,r,col,t,:,:), y.flash_count(n,r,col,t,:,:)] ...
						= calculate_mean_spiking(spike_count_current, flash_pos_current, time_windows(t,:));
					[y.mean_spiking_ds(n,r,col,t,:,:), y.flash_count_ds(n,r,col,t,:,:)] ...
						= calculate_mean_spiking(spike_count_current_ds, flash_pos_ds_current, time_windows(t,:));

					% Continue with analysis if criterion is met
					if sum(spike_count_current_ds) > 0
						
						%%% Grow receptive fields using mutual information
						[y.RF_mask(n,r,col,t,1,:,:), y.RF_ROC(n,r,col,t,1), y.RF_MI(n,r,col,t,1)] = ...
							grow_RF(spike_count_current_ds, flash_pos_ds_current, false, window_size);
						[y.RF_mask(n,r,col,t,2,:,:), y.RF_ROC(n,r,col,t,2), y.RF_MI(n,r,col,t,2)] = ...
							grow_RF(spike_count_current_ds, flash_pos_ds_current, true, window_size);
						
						%%% Fit receptive fields using mutual information
						% z_smooth_current = squeeze(sum(z_smooth(:,:,time_windows(t,1):time_windows(t,end)),3));
						% [y.RF_fit_params(n,r,col,t,:), y.RF_gaussian_fit(n,r,col,t,:,:), ~, ~, ~, ~, y.RF_R2(n,r,col,t)] = ...
						% 	fmgaussfit_ball_strikes(x_pos,y_pos,z_smooth_current,flash_count_smooth,options);
						
						%%% Generate a binary mask based on spiking rate percentiles
						y.perc_mask(n,r,col,t,:,:,:) = perc_mask(spike_count_current_ds, flash_pos_ds_current, percentiles);
						% y.smoothed_perc_mask(n,r,col,t,:,:,:) = perc_mask(smoothed_spike_count_current_ds, [], percentiles);
							
					end
				end
			end
			
			%%% Calculate ROC values between different flash types
			% z = [];
			% z.spike_resp = spike_resp_ds;
			% z.flash_pos = flash_pos_ds;
			% z.flash_col = flash_col_ds;
			% z.flash_zone = flash_zone_ds;
			% z.target_flash = target_flash_ds;
			% [y.flash_roc(n,r,:,:), y.comparisons] = calculate_roc(z, time_windows, roc_descriptions);

		end
	end

	% Post-processing
	% Collapse across repetitions for certain outputs before returning
	y.mean_spiking_ds = squeeze(mean(y.mean_spiking_ds, 2));
	y.mean_spiking = squeeze(mean(y.mean_spiking, 2));
	y.flash_count_ds = squeeze(mean(y.flash_count_ds, 2));
	y.flash_count = squeeze(mean(y.flash_count, 2));
	% y.smoothed_perc_mask = squeeze(mean(y.smoothed_perc_mask, 2));
	y.perc_mask = squeeze(mean(y.perc_mask, 2));
	y.RF_mask = squeeze(mean(y.RF_mask, 2));
	y.RF_ROC = squeeze(mean(y.RF_ROC, 2));
	y.RF_MI = squeeze(mean(y.RF_MI, 2));

end


function [mean_spiking, flash_count] = calculate_mean_spiking(spike_resp, flash_pos, time_window)

	% time_window is in milliseconds, converting to seconds for dt
	% for mean spiking to be in Hz
	dt = (1 + time_window(:,2) - time_window(:,1))/1e3;
	mean_spiking = zeros(13,26,'single');
	flash_count = zeros(13,26,'single');

	for x_pos = 1:13
		for y_pos = 1:26

			% Test for bounding boxes
			correct_x = flash_pos(:,1) == x_pos;
			correct_y = flash_pos(:,2) == y_pos;
			spike_pos_inds = find(correct_x & correct_y);

			% Take mean spike rate for that location
			if ~isempty(spike_pos_inds)
				mean_spiking(x_pos,y_pos) = mean(spike_resp(spike_pos_inds)) / dt;
				flash_count(x_pos,y_pos) = length(spike_pos_inds);
			end
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mean_resp, flash_count] = convert_single_trial_to_mean_resp(spike_resp, flash_pos)

	m = 13;
	n = 26;
	trial_length = size(spike_resp, 2);
	mean_resp = zeros(n, m, trial_length);
	flash_count = zeros(n, m);
	for x = 1:m
		for y = 1:n
			ind = find(flash_pos(:,1)==x & flash_pos(:,2)==y);
            if ~isempty(ind)
                mean_resp(y,x,:) = mean(spike_resp(ind,:), 1);
                flash_count(y,x) = length(ind);
            end
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binary_masks = perc_mask(spike_count, flash_pos, percentiles)

	% Get number of percentiles
	num_pcts = length(percentiles);

	% Make data arrays
	spikes_by_position = zeros(13,26,'single');
	binary_masks = zeros(num_pcts,13,26,'logical');

	if all(size(spike_count) > 1)

		% If already in positional form, just calculate the percentiles
		% and store accordingly
		pcts = prctile(spike_count, percentiles, 'all');
		for p = 1:num_pcts
			binary_masks(p,:,:) = spike_count > pcts(p);
		end
		
	else

		% If as a list of spikes, iterate over all locations
		for x_pos = 1:13
			for y_pos = 1:26
				
				% Test for bounding boxes
				correct_x = flash_pos(:,1) == x_pos;
				correct_y = flash_pos(:,2) == y_pos;
				spike_pos_inds = find(correct_x & correct_y);
				
				% In the appropriate locations, add the mean spike count
				spikes_by_position(x_pos,y_pos) = mean(spike_count(:,spike_pos_inds));
			end
		end

		% Take percentiles of the mean spike counts, and store accordingly
		pcts = prctile(spikes_by_position, percentiles, 'all');
		for p = 1:num_pcts
			binary_masks(p,:,:) = spikes_by_position > pcts(p);
		end
	end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RF_mask, ROC, MI] = grow_RF(s, flash_pos, force_positive, window_size)

	% flash_pos should be [1:13] for x, and [1:26] for y
	max_x = double(max(flash_pos(:,1)));
	max_y = double(max(flash_pos(:,2)));

	% intial windows sizes
	if isempty(window_size)
		window_x = 3;
		window_y = 3;
	else
		window_x = window_size(1);
		window_y = window_size(2);
	end
	num_flashes = size(flash_pos, 1);

	% Step 1: find window_x X window_y RF that gives best ROC

	best_ROC = 0;
	best_MI = 0;
	best_x = [];
	best_y = [];
	for x_pos = 1:max_x-window_x+1
		for y_pos = 1:max_y-window_y+1
			ind1 = logical(flash_pos(:,1)>=x_pos & flash_pos(:,1)<x_pos+window_x & flash_pos(:,2)>=y_pos & flash_pos(:,2)<y_pos+window_y);
			ind2 = logical(ones(num_flashes,1) - ind1);
			%         current_ROC = roc_masse(s(ind2)',s(ind1)');
			mask_fraction = window_x*window_y/(max_x*max_y);
			%         current_ROC = threshold_discrimation(s(ind2)',s(ind1)', 0.5, rectify_ROC);
			current_MI = calculate_mutual_info(s(ind2),s(ind1),mask_fraction);
			if force_positive
				if mean(s(ind1)) < mean(s(ind2))
					current_MI = -1;
				end
			end
			%         if rectify_ROC
			%             current_ROC = abs(0.5  - current_ROC);
			%         end
			%         current_MI = convert_ROC_to_MI(current_ROC, mask_fraction);
			if current_MI > best_MI
				%             best_ROC = current_ROC;
				best_MI = current_MI;
				best_x = x_pos;
				best_y = y_pos;
			end
		end
	end

	% STEP 2: Find region with best ROC (non-rectified)
	if ~isempty(best_x)
		current_mask = zeros(max_x, max_y);
		current_mask(best_x:best_x+window_x-1,best_y:best_y+window_y-1) = 1;
		ind_target = logical(flash_pos(:,1)>=best_x & flash_pos(:,1)<best_x+window_x & flash_pos(:,2)>=best_y & flash_pos(:,2)<best_y+window_y);
		continue_search = true;
		
		while continue_search
			
			dx = current_mask(2:end,:) - current_mask(1:end-1,:);
			dx1 = cat(1, zeros(1, max_y), dx);
			dx2 = cat(1, dx, zeros(1, max_y));
			dy = current_mask(:, 2:end) - current_mask(:,1:end-1);
			dy1 = cat(2, zeros(max_x, 1), dy);
			dy2 = cat(2, dy, zeros(max_x, 1));
			search_region = logical(max(0, min(1,abs(dx1) + abs(dx2) + abs(dy1) + abs(dy2)) - current_mask));
			%     figure
			%     subplot(1,2,1),
			%     imagesc(current_mask)
			%     subplot(1,2,2),
			%     imagesc(search_region)
			search_region = find(search_region);
			current_ROC = zeros(1,length(search_region));
			current_MI = zeros(1,length(search_region));
			%     current_MI = zeros(1,length(search_region));
			for j = 1:length(search_region)
				[current_x, current_y] = ind2sub([max_x max_y], search_region(j));
				ind_current = flash_pos(:,1) == current_x & flash_pos(:,2) == current_y;
				ind_new_target = ind_target | ind_current;
				ind_rest = logical(ones(num_flashes,1) - ind_new_target);
				%         current_ROC(j) = roc_masse(s(ind_rest)',s(ind_new_target)');
				mask_fraction = (sum(sum(current_mask))+1)/numel(current_mask);
				%         current_ROC(j) = threshold_discrimation(s(ind_rest)',s(ind_new_target)', 0.5, rectify_ROC);
				%         if rectify_ROC
				%             current_ROC(j) = abs(0.5 - current_ROC(j));
				%         end
				%         current_MI(j) = convert_ROC_to_MI(current_ROC(j), mask_fraction);
				current_MI(j) = calculate_mutual_info(s(ind_rest),s(ind_new_target), mask_fraction);
				if current_MI(j) > best_MI
					break
				end
				%         mask_fraction = (sum(sum(current_mask))+1)/numel(current_mask);
				%         current_MI(j) = convert_ROC_to_MI(current_ROC(j), mask_fraction);
			end
			
			%     [max_ROC, max_ind] = max(current_ROC);
			[max_MI, max_ind] = max(current_MI);
			if max_MI >= best_MI
				%         best_ROC = max_ROC;
				best_ROC = current_ROC(max_ind);
				best_MI = max_MI;
				% add new flash location to mask
				[current_x, current_y] = ind2sub([max_x max_y], search_region(max_ind));
				current_mask(current_x, current_y) = 1;
				% add new trials to target list
				ind_new = flash_pos(:,1) == current_x & flash_pos(:,2) == current_y;
				ind_target = ind_target | ind_new;
			else
				continue_search = false;
				ind_rest = logical(ones(num_flashes,1) - ind_target);
				best_ROC = roc_masse(s(ind_rest)',s(ind_new_target)');
			end
		end
		
		RF_mask = current_mask;
		ROC = best_ROC;
		MI = best_MI;
	else
		MI = NaN;
		ROC = NaN;
		RF_mask = 0;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MI = calculate_mutual_info(s1, s2, p2)

	% assuming that s1 and s2 are INTEGERS,  spike counts for segements 1 and
	% 2, respectively
	% p2 is the probability that flash occured in segement 2

	% MI = sum p(x,y)*log(p(x,y)/p(x)*p(y))
	% where x is grid segment (1 or 2) where the flash occured
	% and y is the the grid segment that we decoded, based on the ROC score

	% multiplying p2 basically chages the effective size of the region we're
	% interested in
	p2 = 2*double(p2);

	n1 = length(s1);
	n2 = length(s2);
	p1 = 1 - p2; % probability that flash occured in segement 1


	m = max(s2);
	MI = zeros(1,m);
	for threshold = 1:m
		% values equal to or above threshold will be decoded as segment 2
		n21 = sum(s1>=threshold);
		n22 = sum(s2>=threshold);
		pd1a1 = p1*(1-n21/n1);
		pd2a1 = p1*n21/n1;
		pd1a2 = p2*(1-n22/n2);
		pd2a2 = p2*n22/n2;
		pd1 = pd1a1 + pd1a2;
		pd2 = pd2a1 + pd2a2;
		% MI = flash occured in 1, decoded in 1 + flash occured in 1, decoded in 2
		% ... + flash occured in 2, decoded in 1 + flash occured in 2, decoded in 2
		MI(threshold) = pd1a1*log2(pd1a1/(p1*pd1)) +  pd2a1*log2(pd2a1/(p1*pd2)) ...
			+  pd1a2*log2(pd1a2/(p2*pd1)) +  pd2a2*log2(pd2a2/(p2*pd2));
	end

	MI = max(MI);
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resp, flash_count_smooth] = smooth_data(mean_spike_count, xx, yy, flash_count)

	f = [[1/(1+sqrt(2)), 1/2, 1/(1+sqrt(2))]; [1/2, 1, 1/2]; [1/(1+sqrt(2)), 1/2, 1/(1+sqrt(2))]];

	mean_spikes_smooth = convn(mean_spike_count, reshape(f, [3,3,1]), 'same');
	flash_count_smooth = conv2(flash_count, f, 'same');

	resp = mean_spikes_smooth ./ flash_count_smooth;

	% Mask the edges of the convolved data to prevent anomalous percentile judgements
	resp(:,1,:) = 0.;
	resp(1,:,:) = 0.;
	resp(:,end,:) = 0.;
	resp(end,:,:) = 0.;

	% % zz = mean_spike_count
	% % zz_smooth = resp
	% [m,n,trial_length] = size(zz);
	% zz_smooth = zeros(m,n,trial_length);
	% flash_count_smooth = zeros(m,n);

	% x = reshape(xx,m*n,1);
	% y = reshape(yy,m*n,1);
	% z = reshape(zz,m*n,trial_length);
	% fc = reshape(flash_count,numel(flash_count),1);
	% xy = [x y];
	% K = 10;

	% for i = 1:m
	% 	for j = 1:n
	% 		current_point = [xx(i,j) yy(i,j)];
	% 		[idx, d] = knnsearch(xy,current_point,'K',K);
	% 		ind = find(d<2);
	% 		total_weight = 0;
	% 		for k = ind
	% 			weight = 1/(1+d(k));
	% 			if ~isnan(z(idx(k)))
	% 				zz_smooth(i,j,:) = zz_smooth(i,j,:) + reshape(z(idx(k),:)*fc(idx(k))*weight,[1 1 trial_length]);
	% 				flash_count_smooth(i,j) = flash_count_smooth(i,j) + fc(idx(k))*weight;
	% 				total_weight = total_weight + fc(idx(k))*weight;
	% 			end
	% 		end
	% 		zz_smooth(i,j,:) = zz_smooth(i,j,:)/total_weight;
	% 	end
	% end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flash_roc, comparisons] = calculate_roc(y, time_windows, description)

	num_neurons = length(y);
	num_windows = size(time_windows, 1);
	flash_roc = NaN*ones(399, num_windows, 'single');
	comparisons = cell(399, 2);

	flash_colors = {'Red';'Green';'White'};

	ind_group = cell(1,24);
	for n = 1:num_neurons
		ind_group{1} = find(y.target_flash == 1);
		ind_group{2} = find(y.flash_zone == 1 & y.flash_col==1);
		ind_group{3} = find(y.flash_zone == 2 & y.flash_col==1);
		ind_group{4} = find(y.flash_zone == 1 & y.flash_col==2);
		ind_group{5} = find(y.flash_zone == 2 & y.flash_col==2);
		ind_group{6} = find(y.flash_zone == 1 & y.flash_col==3);
		ind_group{7} = find(y.flash_zone == 2 & y.flash_col==3);
		ind_group{8} = find(y.target_flash == 0);
		ind_group{9} = find(y.target_flash == 0 & y.flash_col<=2);
		ind_group{10} = find(y.target_flash == 0 & y.flash_col<=2 & y.flash_zone<=4);
		ind_group{11} = find(y.target_flash == 0 & y.flash_col<=2 & y.flash_zone<=2);
		ind_group{12} = find(y.target_flash == 0 & y.flash_col==1);
		ind_group{13} = find(y.target_flash == 0 & y.flash_col==1 & y.flash_zone<=4);
		ind_group{14} = find(y.target_flash == 0 & y.flash_col==1 & y.flash_zone<=2);
		ind_group{15} = find(y.target_flash == 0 & y.flash_col==2);
		ind_group{16} = find(y.target_flash == 0 & y.flash_col==2 & y.flash_zone<=4);
		ind_group{17} = find(y.target_flash == 0 & y.flash_col==2 & y.flash_zone<=2);
		ind_group{18} = find(y.target_flash == 0 & y.flash_col==3);
		ind_group{19} = find(y.target_flash == 0 & y.flash_col==3 & y.flash_zone<=4);
		ind_group{20} = find(y.target_flash == 0 & y.flash_col==3 & y.flash_zone<=2);
		ind_group{21} = find(y.flash_col == 3 & (y.flash_zone==1 | y.flash_zone==4));
		ind_group{22} = find(y.flash_col == 3 & (y.flash_zone==2 | y.flash_zone==3));
		ind_group{23} = find(y.flash_col == 3 & (y.flash_zone==2 | y.flash_zone==3 | y.flash_zone==4));
		ind_group{24} = find(y.flash_col == 3 & (y.flash_zone==1 | y.flash_zone==3 | y.flash_zone==4));
	end

	for t = 1:num_windows
		comparison_count = 0;
		for col1 = 1:3
			for zone1 = 1:5
				for col2 = 1:3
					for zone2 = 1:5
						if ~(col1==col2 & zone1==zone2)
							comparison_count = comparison_count+1;
							ind1 = find(y.flash_col == col1 & y.flash_zone == zone1);
							ind2 = find(y.flash_col == col2 & y.flash_zone == zone2);
							comparisons{comparison_count, 1} = [flash_colors{col1} num2str(zone1)];
							comparisons{comparison_count, 2} = [flash_colors{col2} num2str(zone2)];
							if ~isempty(ind1) & ~isempty(ind2)
								s1 = sum(y.spike_resp(ind1,time_windows(t,1):time_windows(t,end)),2);
								s2 = sum(y.spike_resp(ind2,time_windows(t,1):time_windows(t,end)),2);
								flash_roc(comparison_count,t) = roc_masse(s1,s2);
							end
						end
					end
				end
			end
		end
		
		% calculate roc for combined groups
		for i = 1:7
			for j = 4:24
				comparison_count = comparison_count+1;
				comparisons{comparison_count, 1} = description{i};
				comparisons{comparison_count, 2} = description{j};
				if ~isempty(ind_group{n,i}) & ~isempty(ind_group{n,j})
					s1 = sum(y.spike_resp(ind_group{n,i},time_windows(t,1):time_windows(t,end)),2);
					s2 = sum(y.spike_resp(ind_group{n,j},time_windows(t,1):time_windows(t,end)),2);
					flash_roc(comparison_count,t) = roc_masse(s1,s2);
				end
			end
		end
	end


	comparisons = comparisons(1:comparison_count, :);
	flash_roc = flash_roc(1:comparison_count,:);
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mem_sac_response, mem_sac_pref_angle, mem_sac_p_val, mem_sac_mod_depth, epochs]  = get_mem_saccade_response(session)


	epochs = [-499 0;50 100; 75 125; 75 150; 50 150; 100 200; 100 300; 500 1000; 500 1300; 1300 1600] + 500;

	sorted_window = get_sorted_window(session.x);
	num_neurons = length(session.x(1).spike_times);
	mem_sac_response = zeros(num_neurons, 8, 2200);
	count = zeros(num_neurons,8);

	num_epochs = size(epochs, 1);
	mem_sac_pref_angle = zeros(num_neurons, num_epochs);
	mem_sac_p_val = ones(num_neurons, num_epochs);
	mem_sac_mod_depth = zeros(num_neurons, num_epochs);

	for n = 1:num_neurons
		sac_direction = [];
		spike_count = [];
		for i = 1:length(session.x)
			if strcmp(session.x(i).trial_type,'memory_saccade') & ~isempty(session.x(i).reward_time)
				cond = session.x(i).condition_number;
				
				if i>=sorted_window(n,1) && i<=sorted_window(n,2)
					sac_direction = [sac_direction; cosd((cond-1)*45) sind((cond-1)*45)];
					spike_count_current = zeros(1, num_epochs);
					st = session.x(i).spike_times{n} - session.x(i).sac_target_on + 500;
					for k = 1:num_epochs
						spike_count_current(k) = sum(st>epochs(k, 1) & st<=epochs(k, 2));
					end
					spike_count = [spike_count; spike_count_current];
					count(n,cond) = count(n,cond) + 1;
					st = st(st>0 & st<= 2200);
					mem_sac_response(n,cond,st) = mem_sac_response(n,cond,st) + 1;
				end
			end
		end
		if size(spike_count, 1) >= 40 % need at least 40 trials
			for k = 1:num_epochs
				b = glmfit(sac_direction, spike_count(:,k));
				theta = atan2(b(3),b(2));
				mem_sac_pref_angle(n, k) = theta;
				r = [cos(theta) sin(theta); -sin(theta) cos(theta)];
				sac_direction_rotated = sac_direction*r';
				[b,~,stats] = glmfit(sac_direction_rotated, spike_count(:,k));
				mem_sac_p_val(n, k) = stats.p(2);
				z = [ones(size(spike_count,1),1) sac_direction_rotated];
				H = z\spike_count(:,k);
				spike_count_predicted = z*H;
				spike_err = spike_count_predicted - spike_count(:,k);
				Q = spike_err'*spike_err/size(z,1);
				mem_sac_mod_depth(n, k) = sqrt(sum(H(2:3).^2))/sqrt(Q);
			end
		end
	end
	mem_sac_response = single(mem_sac_response./repmat(count,[1 1 2200]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [green_dp, red_dp, time_windows] = calculate_detect_prob(session, desired_block, time_windows, green_or_cyan)

	if isempty(time_windows)
		time_windows = [[-300:5:300]' [-200:5:400]'];
	end
	sorted_window = get_sorted_window(session.x);
	num_windows = length(time_windows);
	num_trials = length(session.x);
	num_neurons = length(session.x(1).spike_times);
	green_dp = zeros(num_neurons,num_windows);
	red_dp = zeros(num_neurons, num_windows);

	for t = 1:num_windows
		for k = 1:num_neurons
			green_correct = [];
			red_correct = [];
			green_failed = [];
			red_failed = [];
			for i = 1:num_trials
				if i>=sorted_window(k,1) && i<=sorted_window(k,2)
					if strcmp(session.x(i).trial_type,'ball_strikes') & session.x(i).block_number == desired_block & ~isempty(session.x(i).target_flash) ...
							& session.x(i).TrialError <= 1
						
						% just confirming whether it's a red or green target
						% flahs based on its position
						if ~isempty(session.x(i).green_target_pos)
							flash_pos = session.x(i).green_target_pos;
						elseif ~isempty(session.x(i).red_target_pos)
							flash_pos = session.x(i).red_target_pos;
						end
						target_zone = get_zone_precise(flash_pos(end,:),green_or_cyan);
						if target_zone > 2
							error('Something went wrong')
						end
						if (target_zone == 1 & ~isempty(session.x(i).green_target_pos)) | (target_zone == 2 & ~isempty(session.x(i).red_target_pos))
							error('Something went wrong')
						end
						
						st = session.x(i).spike_times{k} - session.x(i).target_flash;
						if target_zone==2 & session.x(i).TrialError == 0
							green_correct = [green_correct; sum(st>=time_windows(t,1) & st<=time_windows(t,2))];
						elseif target_zone==1 & session.x(i).TrialError == 0
							red_correct = [red_correct; sum(st>=time_windows(t,1) & st<=time_windows(t,2))];
						elseif target_zone==2 & session.x(i).TrialError == 1
							green_failed = [green_failed; sum(st>=time_windows(t,1) & st<=time_windows(t,2))];
						elseif target_zone==1 & session.x(i).TrialError == 1
							red_failed = [red_failed; sum(st>=time_windows(t,1) & st<=time_windows(t,2))];
						end
					end
				end
			end
			
			green_dp(k, t) = roc_masse(green_failed, green_correct);
			red_dp(k, t) = roc_masse(red_failed, red_correct);
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_resp_ds, flash_col_ds, flash_pos_ds, flash_zone_ds, target_flash_ds] = downsample_target_flashes_combined(spike_resp, flash_pos, flash_col, trial_type)

	spike_resp_ds = [];
	flash_pos_ds = [];
	flash_col_ds = [];
	target_flash_ds = [];
	flash_zone_ds = [];
	falsh_count_hist = [];
	if ~isempty(flash_pos)
		% calculate the flash count histogram
		for x = 1:13
			for y = 1:26
				% count flash map flash counts (trial type = 1)
				ind = find(flash_col == 3 & flash_pos(:,1)==x & flash_pos(:,2)==y & trial_type==1 & flash_col == 3);
				falsh_count_hist = cat(1, falsh_count_hist, length(ind));
			end
		end
		
		
		
		for tt = 1:2
			% select only a subset of  flahses based on the flash count histogram
			% calculated for non-target flahses
			for x = 1:13
				for y = 1:26
					
					q = randperm(length(falsh_count_hist));
					flahses_to_use = falsh_count_hist(q(1));
					ind = find(flash_pos(:,1)==x & flash_pos(:,2)==y & trial_type==tt & flash_col == 3);
					if length(ind) >= flahses_to_use
						q = randperm(length(ind));
						ind = ind(q(1:flahses_to_use));
					else
						% if there are not enough flahses, use flashes from
						% adjacent points
						num_additional_flashes = flahses_to_use - length(ind);
						ind2 = find(flash_pos(:,1)>=x-1 & flash_pos(:,1)<=x+1 & flash_pos(:,2)>=y-1 ...
							& flash_pos(:,2)<=y+1 & trial_type==tt & flash_col == 3);
						ind2 = setdiff(ind2, ind);
						if length(ind2)>num_additional_flashes
							q = randperm(length(ind2));
							ind2 = ind2(q(1:num_additional_flashes));
						end
						ind = union(ind, ind2);
					end
					flash_pos_ds = cat(1, flash_pos_ds, repmat([x y],length(ind),1));
					% hijacking color, replacing with trial type!!!
					flash_col_ds = cat(1, flash_col_ds, repmat(tt,length(ind),1));
					target_flash_ds = cat(1, target_flash_ds, ones(length(ind),1));
					flash_zone_ds = cat(1, flash_zone_ds, repmat(3,length(ind),1));
					spike_resp_ds = cat(1, spike_resp_ds, spike_resp(ind,:));
					
				end
			end
		end
	end
end
