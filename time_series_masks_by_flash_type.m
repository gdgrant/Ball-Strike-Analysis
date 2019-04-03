function time_series_masks_by_flash_type()

	% load('./result_data/result_data_old_code.mat')
	% load('./result_data/metadata_old_code.mat')

	% num_neurons = size(results.metadata.iso_qual, 1);
	% fn = fieldnames(results.metadata);
	% for i=1:length(fn)

	% 	if any(size(results.metadata.(fn{i})) == 0) && ~any(size(metadata.(fn{i})) == 0)
	% 		size(results.metadata.(fn{i}))
	% 		size(metadata.(fn{i}))
	% 		results.metadata.(fn{i}) = metadata.(fn{i})(1:num_neurons,:);
	% 	end
	% end

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

	names = {'All', 'PPC', 'PFC (all)', 'PFC (wave>=10)', 'PFC (wave<10)'};
	sets  = {all_neurons, PPC_neurons, PFC_neurons, PFC_neurons_wave_g215, PFC_neurons_wave_l215};
	
	subject = 'Spiking';

	for i=2:3

		% results.data.perc_mask : neurons x colors x windows x percs x x_locs x y_locs

		switch subject
			case 'Perc'
				collection = results.data.perc_mask(sets{i},:,:,perc,:,:);
				mask_name = ['Perc Mask (' num2str(percentiles(perc)) '%)'];
			case 'SPerc'
				collection = results.data.smoothed_percentiles(sets{i},:,:,perc,:,:);
				mask_name = ['Smoothed Perc Mask (' num2str(percentiles(perc)) '%)'];
			case 'RF'
				% plc = squeeze(mean(results.data.ROCs.RF_mask, 2));
				% collection = plc(sets{i},:,:,1,:,:);
				collection = results.data.ROCs.RF_mask(sets{i},:,:,1,:,:);
				mask_name = 'RF Mask (fp=False)';
			case 'Spiking'
				collection = results.data.mean_spiking(sets{i},:,:,:,:);
				% collection_count = results.data.flash_count(sets{i},:,:,:,:);
				mask_name = 'Mean Spiking';
			case 'Spiking DS'
				collection = results.data.mean_spiking_ds(sets{i},:,:,:,:);
				mask_name = 'Mean Spiking DS';
		end

		% Remove dimensions specific to data case
		collection = squeeze(collection);

		% ind to pos : x=xind-1, y=yind-13.5
		% pos to ind : xind=x+1, yind=y+13.5
		% red_zone   : x=[4 8],  y=[0.5 4.5]   --> [5,9] [14,18]
		% green_zone : x=[4 8],  y=[-4.5 -0.5] --> [5,9] [9,13]

		x_inds = 1:13;
		y_inds = 1:26;
		[x_grid, y_grid] = meshgrid(x_inds, y_inds);

		red_zone   = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 14 & y_grid <= 18);
		green_zone = (x_grid >= 5 & x_grid <= 9) & (y_grid >= 9 & y_grid <= 13);
		white_zone = ~(red_zone | green_zone);

		red_flashes   = squeeze(collection(:,1,:,:,:));
		green_flashes = squeeze(collection(:,2,:,:,:));
		white_flashes = squeeze(collection(:,3,:,:,:));

		colors = {'Red', 'Green', 'White'};
		flashes = {red_flashes, green_flashes, white_flashes};
		zones = {red_zone', green_zone', white_zone'};

		ax = []; 
		figure;
		for c=1:3

			% neurons x windows x locs
			red_zone_response_by_neuron   = flashes{c}(:,:,zones{1});
			green_zone_response_by_neuron = flashes{c}(:,:,zones{2});
			white_zone_response_by_neuron = flashes{c}(:,:,zones{3});

			% neurons x windows
			amb = squeeze(mean(red_zone_response_by_neuron-green_zone_response_by_neuron, 3));
			[h, p] = ttest(amb);

			% windows
			red_zone_traj   = squeeze(mean(red_zone_response_by_neuron, [1,3]));
			green_zone_traj = squeeze(mean(green_zone_response_by_neuron, [1,3]));
			white_zone_traj = squeeze(mean(white_zone_response_by_neuron, [1,3]));

			windows = 1:31;

			ax1 = subplot(1,3,c);
			plot(...
				windows, red_zone_traj,   'r', ...
				windows, green_zone_traj, 'g', ...
				windows, white_zone_traj, 'k');
			hold on

			%legend({'Red Zone';'Green Zone';'White Zone'});

			for w=windows
				if p(w)<0.1
					plot([w,w+1], 0.09*ones(1,2), 'k-')
				end
			end
			hold off

			ax = cat(1,  ax, ax1);
			set(gca, 'xlim', [1,31]);
			xlabel('Time Window');
			ylabel('Mean Mask Value');
			title(strcat(colors{c}, ' Flash Resp'));
		end
		linkaxes(ax, 'y');
		% sgtitle(['Time Series Mean Mask Value Per Flash Type, ' names{i} ' Neurons']);
		sgtitle([mask_name ' Value Per Flash Type, ' names{i} ' Neurons']);
	end
end