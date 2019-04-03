function [valid_flash_inds] = cull_flashes(trial, num_dead_flashes)
% Use flashes up to target flash, or second to last flash

	% Set num_dead_flashes to a default value
	if nargin < 2
		num_dead_flashes = 1;
	end

	% Test for a target flash
	zone1 = trial.flash_color==1 & trial.flash_zone==1;
	zone2 = trial.flash_color==2 & trial.flash_zone==2;
	target_flash_ind = find(zone1 | zone2);

	% If no target flash, use second to last flash
	if isempty(target_flash_ind)
		last_flash = length(trial.flash_time) - 1;
	else
		last_flash = target_flash_ind;
	end

	% Generate valid flash indices
	valid_flash_inds = num_dead_flashes+1:last_flash;
	
end

