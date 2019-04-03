function [spike_resp_ds, flash_col_ds, flash_pos_ds, flash_zone_ds, target_flash_ds] = downsample_target_flashes(spike_resp, flash_pos, flash_col, flash_zone, green_or_cyan, sample_white)

	%%% Use all the flashes not in the target zone to develop a flash histogram
	%%% Input is a single neuron for an entire session, plus parameters

	% Set up data collection
	spike_resp_ds = [];
	flash_pos_ds = [];
	flash_col_ds = [];
	target_flash_ds = [];
	flash_zone_ds = [];

	% Designate default sampling method
	if ~exist('sample_white', 'var') || isempty(sample_white)
		sample_white = false;
	end

	% Exit the function if there are no flashes
	if isempty(flash_pos)
		return
	end

	% Iterate over the three colors.  At the end of each pass,
	% add the desired flashes to the _ds arrays


	if ~sample_white
		col = 3;
		
		ind = find(flash_col == col);
		num_flashes = length(ind);
		
		flash_pos_ds = cat(1, flash_pos_ds, [flash_pos(ind,1) flash_pos(ind,2)]);
		flash_col_ds = cat(1, flash_col_ds, repmat(col, num_flashes, 1));
		target_flash_ds = cat(1, target_flash_ds, ones(num_flashes, 1));
		spike_resp_ds = cat(1, spike_resp_ds, spike_resp(ind,:));
		flash_zone_ds = cat(1, flash_zone_ds, flash_zone(ind,:));
	end


	if ~sample_white
		max_col = 2;
	else
		max_col = 3;
	end

	for col = 1:max_col
		
		% Start the flash count histogram
		% This histogram records the number of flashes in each location and condition
		flash_count_hist = [];
		
		% Build histogram
		% Iterate over all positions
		for x = 1:13
			for y = 1:26
				
				% Identify the current position's zone : 1 (Red), 2 (Green), or 3 (White)
				current_zone = get_zone(x-1, y-13.5, green_or_cyan);
				
				% If the current zone doesn't match the desired color, or if the desired
				% color is white, append the tally of relevant flashes to the histogram
				if current_zone ~= col || col == 3
					
					% Obtain flashes relevant to this position and color
					ind = find(flash_col == col & flash_pos(:,1)== x & flash_pos(:,2) == y);
					
					% Append flash tally
					flash_count_hist = cat(1, flash_count_hist, length(ind));
				end
			end
		end
		
		% If opting to not sample white, and if the color is white,
		% skip the sampling step (white flashes are already uniformly distributed)
		if ~sample_white && col==3
			continue
		end
		
		% Sample from histogram
		% Iterate over all positions
		ind_list = [];	% Keep track of inds over positions
		for x = randperm(13)
			for y = randperm(26)
				
				% Identify the current position's zone : 1 (Red), 2 (Green), or 3 (White)
				current_zone = get_zone(x-1, y-13.5, green_or_cyan);
				
				% Select a subset of flashes based on the flash count histogram, as
				% calculated for non-target flashes
				flashes_to_use = flash_count_hist(randi(length(flash_count_hist)));
				
				% Obtain flashes relevant to this position and color
				% that have not yet been used
				ind = find(flash_col == col & flash_pos(:,1)==x & flash_pos(:,2)==y);
				ind = setdiff(ind, ind_list);
				
				% If the number of available flashes is greater than the number of flashes
				% requested to do the sampling, simply pull those flashes.  Otherwise, use
				% flashes from adjacent points to supplement.
				if length(ind) >= flashes_to_use
					
					% Take 'flashes_to_use' flash samples from the index array
					ind = ind(randperm(length(ind), flashes_to_use));
					
				else
					
					% Determine the number of extra flashes required to supplement
					num_additional_flashes = flashes_to_use - length(ind);
					
					% Obtain those flashes and remove already-used flashes
					ind2 = find(flash_col == col & flash_zone == current_zone ...
						& flash_pos(:,1)>=x-1 & flash_pos(:,1)<=x+1 ...
						& flash_pos(:,2)>=y-1 & flash_pos(:,2)<=y+1);
					
					% Cull to avoid repetition and reduce number of samples to required number
					ind2 = setdiff(ind2, union(ind, ind_list));
					
					% Cap the number of extra flash indices to the number of required
					% extra flashes, then sample among those indices
					num_inds = min(length(ind2), num_additional_flashes);
					ind2 = ind2(randperm(num_inds));
					
					% Combine the extra flashes with the already-relevant flashes
					ind = union(ind, ind2);
				end

				if ~iscolumn(ind)
					ind = ind';
				end
				ind_list = cat(1, ind_list, ind);
				
				% Append flash data to the appropriate arrays
				num_flashes = length(ind);
				flash_pos_ds = cat(1, flash_pos_ds, repmat([x y],num_flashes,1));
				flash_col_ds = cat(1, flash_col_ds, repmat(col,num_flashes,1));
				target_flash_ds = cat(1, target_flash_ds, ones(num_flashes,1));
				flash_zone_ds = cat(1, flash_zone_ds, repmat(current_zone,num_flashes,1));
				
			end
		end
		
		spike_resp_ds = cat(1, spike_resp_ds, spike_resp(ind_list,:));
	end
end


function zone = get_zone(x, y, green_or_cyan)

	% RED Target   - Zone 1
	% GREEN Target - Zone 2
	% OUTSIDE      - Zone 3

	if nargin < 3
		green_or_cyan = 'green';
	end

	red_zone_x = [4 8];
	red_zone_y = [0.5 4.5];

	if strcmp(green_or_cyan,'green')
		green_zone_x = [4 8];
		green_zone_y = [-4.5 -0.5];
	elseif strcmp(green_or_cyan,'cyan')
		green_zone_x = [6 10];
		green_zone_y = [5.5 9.5];
	end

	red_x = (red_zone_x(1) <= x) & (x <= red_zone_x(2));
	red_y = (red_zone_y(1) <= y) & (y <= red_zone_y(2));
	green_x = (green_zone_x(1) <= x) & (x <= green_zone_x(2));
	green_y = (green_zone_y(1) <= y) & (y <= green_zone_y(2));

	if red_x && red_y
		zone = 1;
	elseif green_x && green_y
		zone = 2;
	else
		zone = 3;
	end
end