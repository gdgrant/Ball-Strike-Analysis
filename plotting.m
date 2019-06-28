function plotting()

	% files = dir('./result_data/result_data_with_ROC_test.mat');

	% collection = 0;
	% total_neurons = 0;
	% for i = 1:length(files)
	% 	load(['./result_data/' files(i).name])
	% 	neurons = size(results.data.perc_mask, 1);
		
	% 	% Collapse across neurons when loaded, and take weighted mean
	% 	% to get the right statistics
	% 	collection = collection + neurons*mean(results.data.perc_mask, 1);
	% 	total_neurons = total_neurons + neurons;
	% end

	% % colors x windows x percs x x_locs x y_locs
	% collection = squeeze(collection ./ total_neurons);
	% percentiles = [80,85,90,95];
	% win = 20;

	% name = {'Red' 'Green' 'White'};
	% for p=1:4
	% 	figure
	% 	for i=1:3
	% 		plot_data = squeeze(collection(i,win,p,:,:));
	% 		subplot(1,3,i),imagesc(0:12,-12:5:12.5,plot_data');
	% 		set(gca,'YDir','normal')
	% 		axis equal;
	% 		axis tight;
	% 		colorbar;
	% 		xlabel('X-location')
	% 		ylabel('Y-location')
	% 		title(strcat(name{i}, ' Spike Focus'))
	% 	end
	% 	sgtitle([num2str(percentiles(p)) 'th Percentile, ' num2str(win) 'th Window'])
	% end


% end

	load('./result_data/result_data.mat')
	win = 20;

	% Neuron channels and iso qual conditions
	PPC_region       = results.metadata.channel_num <= 32;
	PFC_region       = results.metadata.channel_num > 32;
	iso_qual_region  = results.metadata.iso_qual >= 3;
	waveform_regionA = results.metadata.waveform_width >= 10;
	waveform_regionB = results.metadata.waveform_width < 10;

	all_neurons = 1:size(results.data.perc_mask,1);
	PPC_neurons = logical(PPC_region .* iso_qual_region);
	PFC_neurons = logical(PFC_region .* iso_qual_region);

	for c=2:3

		switch c
			case 1
				collection = squeeze(mean(results.data.mean_spiking, [1]));
				% collection = squeeze(mean(results.data.mean_spiking(PFC_neurons,:,:,:,:), [1]));
				ttl = 'Mean Spiking (all)';
				plot_data = squeeze(collection(:,win,:,:));
			case 2
				collection = squeeze(mean(results.data.mean_spiking_ds, [1]));
				ttl = 'Mean Spking DS';
				plot_data = squeeze(collection(:,win,:,:));
			case 3
				collection = squeeze(mean(results.data.perc_mask, [1]));
				ttl = 'Percentile Mask';
				plot_data = squeeze(collection(:,win,4,:,:));
			case 4
				collection = squeeze(mean(results.data.smoothed_percentiles, [1]));
				ttl = 'Smoothed Percentile Mask';
				plot_data = squeeze(collection(:,win,4,:,:));
		end

		name = {'Red' 'Green' 'White'};
		max_color = max(plot_data, [], 'all')
		figure
		for i=1:3
			subplot(1,3,i),imagesc(0:12,-12.5:12.5,squeeze(plot_data(i,:,:))');
			set(gca,'YDir','normal');
			axis equal;
			axis tight;
			colorbar;
			% caxis([0,max_color]);
			xlabel('X-location')
			ylabel('Y-location')
			% title(strcat(name{i}, ' Spike Focus'))
			title(name{i})
		end
		sgtitle([ttl ' | ' num2str(win) 'th Window'])
		end
	end