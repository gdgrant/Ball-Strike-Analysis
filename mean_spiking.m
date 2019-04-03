function mean_spiking()

	load('./result_data/result_data.mat')

	num_neurons = size(results.metadata.iso_qual, 1);
	percentiles = [80,85,90,95];
	perc = 4;
	
	% Neuron channels and iso qual conditions
	PPC_region       = results.metadata.channel_num <= 32;
	PFC_region       = results.metadata.channel_num > 32;
	iso_qual_region  = results.metadata.iso_qual >= 3;
	waveform_regionA = results.metadata.waveform_width >= 10;
	waveform_regionB = results.metadata.waveform_width < 10;

	active_region = squeeze(mean(results.data.mean_spiking, [2,3,4,5]))*1e3 > 1.0;
	% active_region = logical(ones(num_neurons,1));

	all_neurons = active_region;
	PPC_neurons = logical(PPC_region .* iso_qual_region .* active_region);
	PFC_neurons = logical(PFC_region .* iso_qual_region .* active_region);
	PFC_neurons_wave_g215 = logical(PFC_neurons .* waveform_regionA .* active_region);
	PFC_neurons_wave_l215 = logical(PFC_neurons .* waveform_regionB .* active_region);

	x_inds = 1:13;
	y_inds = 1:26;
	[x_grid, y_grid] = meshgrid(x_inds, y_inds);

	red_zone   = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 14 & y_grid <= 18);
	green_zone = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 9 & y_grid <= 13);
	white_zone = ~(red_zone | green_zone);

	names = {'All', 'PPC', 'PFC (all)', 'PFC (wave>=10)', 'PFC (wave<10)'};
	sets  = {all_neurons, PPC_neurons, PFC_neurons, PFC_neurons_wave_g215, PFC_neurons_wave_l215};

	colors = {'Red', 'Green', 'White'};
	zones = {red_zone', green_zone', white_zone'};
	windows = 1:31;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

	neuron_set = 3;

	% size: neurons x colors x windows x x_locs x y_locs
	spiking    = results.data.mean_spiking(sets{neuron_set},:,:,:,:);
	spiking_ds = results.data.mean_spiking_ds(sets{neuron_set},:,:,:,:);
	counts     = results.data.flash_count(sets{neuron_set},:,:,:,:);
	counts_ds  = results.data.flash_count_ds(sets{neuron_set},:,:,:,:);

	spiking_by_color    = separate_by_color(spiking);
	spiking_ds_by_color = separate_by_color(spiking_ds);
	counts_by_color     = separate_by_color(counts);
	counts_ds_by_color  = separate_by_color(counts_ds);

	spiking_data = {spiking_by_color;spiking_ds_by_color};
	counts_data = {counts_by_color;counts_ds_by_color};
	name_data = {'(Standard)', '(Downsampled)'};

	for ds=1:2
		figure;
		ax = [];
		for c=1:3

			% Spike responses and flash counts
			% neurons x windows x locs
			red_zone_spk = spiking_data{ds}{c}(:,:,zones{1});
			grn_zone_spk = spiking_data{ds}{c}(:,:,zones{2});
			wht_zone_spk = spiking_data{ds}{c}(:,:,zones{3});

			red_zone_cnt = counts_data{ds}{c}(:,:,zones{1});
			grn_zone_cnt = counts_data{ds}{c}(:,:,zones{2});
			wht_zone_cnt = counts_data{ds}{c}(:,:,zones{3});

			% >>> windows
			red_zone = weighted_func(red_zone_spk, red_zone_cnt);
			grn_zone = weighted_func(grn_zone_spk, grn_zone_cnt);
			wht_zone = weighted_func(wht_zone_spk, wht_zone_cnt);

			ax1 = subplot(1,3,c);
			plot(...
				windows, red_zone, 'r', ...
				windows, grn_zone, 'g', ...
				windows, wht_zone, 'k');
			% hold on

			legend({'Red Zone';'Green Zone';'White Zone'});
			title(strcat(colors{c}, ' Flash Resp'))
			set(gca, 'xlim', [1,31]);
			xlabel('Time Window');
			ylabel('Mean Mask Value');
			ax = cat(1, ax, ax1);
		end
		linkaxes(ax, 'y')
		sgtitle(['Mean Spiking ' name_data{ds} ' Per Flash Type, ' names{neuron_set} ' Neurons']);
	end
end



function y = separate_by_color(x)

	red   = squeeze(x(:,1,:,:,:));
	green = squeeze(x(:,2,:,:,:));
	white = squeeze(x(:,3,:,:,:));

	y = {red; green; white};

end


function z = weighted_func(spike_resp, flash_count)

	% neurons x windows x locs >>> windows
	s = sum(spike_resp .* flash_count, [1,3]) ./ sum(flash_count, [1,3]);
	z = squeeze(s);

end