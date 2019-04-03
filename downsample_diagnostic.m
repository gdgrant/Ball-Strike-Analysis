function downsample_diagnostic

	load('./spatial_category_exp/quincy-13-Mar-2015.mat')
	num_neurons = length(session.x(1).spike_times);
	time_window = [-100 300];

	sorted_window = get_sorted_window(session.x);

	[ball_strikes_response, single_flash_resp, ~, ~, ball_strikes_flash_count] = ...
		calculate_flash_response(session, 1, 1, 0, 'green', time_window, [], [], 'ball_strikes');

	n = 1;
	spike_resp = single_flash_resp(n).spike_resp;
	flash_pos = single_flash_resp(n).flash_pos;
	flash_col = single_flash_resp(n).flash_col;
	flash_zone = single_flash_resp(n).flash_zone;

	reps = 2000;
	count = zeros(reps,3,13,26);
	count_ds = zeros(reps,3,13,26);

	for r=1:reps;
		if mod(r,10) == 0
			disp(['Rep ', num2str(r)])
		end

		[spike_resp_ds, flash_col_ds, flash_pos_ds, flash_zone_ds, target_flash_ds] = ...
			downsample_target_flashes(spike_resp, flash_pos, flash_col, flash_zone, 'green');

		for c=1:3
			for x=1:13
				for y=1:26

					ind = flash_col == c & flash_pos(:,1) == x & flash_pos(:,2) == y;
					count(r,c,x,y) = sum(ind, 'all');

					ind_ds = flash_col_ds == c & flash_pos_ds(:,1) == x & flash_pos_ds(:,2) == y;
					count_ds(r,c,x,y) = sum(ind_ds, 'all');
				
				end
			end
		end
	end

	x_inds = 1:13;
	y_inds = 1:26;
	[x_grid, y_grid] = meshgrid(x_inds, y_inds);

	red_zone   = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 14 & y_grid <= 18);
	% red_zone   = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 9 & y_grid <= 13);
	green_zone = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 9 & y_grid <= 13);
	white_zone = ~(red_zone | green_zone);

	red_flash_target = count(:,1,red_zone);
	red_flash_nontarget = count(:,1,~red_zone);
	red_flash_target_ds = count_ds(:,1,red_zone);
	red_flash_nontarget_ds = count_ds(:,1,~red_zone);

	green_flash_target = count(:,1,green_zone);
	green_flash_nontarget = count(:,1,~green_zone);
	green_flash_target_ds = count_ds(:,1,green_zone);
	green_flash_nontarget_ds = count_ds(:,1,~green_zone);

	ctrs = [0:20];
	[red_target_counts, red_target_centers] = hist(red_flash_target_ds(:),ctrs);
	[red_nontarget_counts, red_nontarget_centers] = hist(red_flash_nontarget_ds(:),ctrs);
	[green_target_counts, green_target_centers] = hist(green_flash_target_ds(:),ctrs);
	[green_nontarget_counts, green_nontarget_centers] = hist(green_flash_nontarget_ds(:),ctrs);

	figure
	plot( ...
		red_target_centers, red_target_counts ./ sum(red_target_counts, 'all'), 'r', ...
		red_nontarget_centers, red_nontarget_counts ./ sum(red_nontarget_counts,'all'), 'b');
	legend({'Red Zone';'Non-Red Zone'});
	xlabel('Number of DSd Flashes Per Location');
	ylabel('Perc. Counts');
	title([num2str(reps) ' Reps'])

	disp(['Red Zone Mean/Var     ' num2str(round(mean(red_flash_target_ds,'all'),3)), '   ', num2str(round(var(red_flash_target_ds, 0, 'all'), 3))])
	disp(['Non-Red Zone Mean/Var ' num2str(round(mean(red_flash_nontarget_ds,'all'), 3)), '   ', num2str(round(var(red_flash_nontarget_ds, 0, 'all'), 3))])


end 