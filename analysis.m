function results = analysis()

	% Get support code and set up environments
	addpath('./support_code/')
	datadir = './spatial_category_exp/';
	savedir = './result_data/';

	% Monkey names
	monkeys = {'quincy';'wahwah'};

	% Make data and metadata fields
	data = [];
	metadata = [];

	data.perc_mask = [];
	data.time_windows = [];
	data.percentiles = [];
	% data.smoothed_percentiles = [];

	data.ROCs = [];
	data.ROCs.RF_mask = [];
	data.ROCs.RF_ROC = [];
	data.ROCs.RF_MI = [];

	data.mean_spiking = [];
	data.mean_spiking_ds = [];
	data.flash_count = [];
	data.flash_count_ds = [];

	metadata.monkey = [];
	metadata.session_date = [];
	metadata.channel_num = [];
	metadata.iso_qual = [];
	metadata.depth = [];
	metadata.waveform_width = [];
	metadata.num_flashes = [];

	% Iterate over monkeys
	for m = 1:length(monkeys);
		
		% Get the desired data associated with the selected monkey
		if strcmp(monkeys{m}, 'wahwah')
			data_files = dir([datadir 'wahwah*Oct*.mat']);
		elseif strcmp(monkeys{m}, 'quincy')
			data_files = dir([datadir 'quincy*.mat']);
		end
		
		% Iterate over the selected files
		for i = 1:length(data_files);
			
			% Load and process a session
			load([datadir data_files(i).name])
			disp(['Session ' data_files(i).name ' loaded'])
			
			% Run analysis on the session
			y = ball_strike_single_unit_stats_single_session_downsample(session);
			neurons = y.num_neurons;

			% Collect data and metadata
			data.perc_mask            = cat(1, data.perc_mask, y.bs_perc_mask1);
			% data.smoothed_percentiles = cat(1, data.smoothed_percentiles, y.bs_smoothed_perc_mask1);
			data.time_windows         = cat(1, data.time_windows, repmat({y.time_windows1}, neurons, 1));
			data.percentiles          = cat(1, data.percentiles, repmat(y.percentiles, neurons, 1));
			data.ROCs.RF_mask         = cat(1, data.ROCs.RF_mask, y.bs_RF_mask1);
			data.ROCs.RF_ROC          = cat(1, data.ROCs.RF_ROC, y.bs_RF_ROC1);
			data.ROCs.RF_MI           = cat(1, data.ROCs.RF_MI, y.bs_RF_MI1);
			data.mean_spiking         = cat(1, data.mean_spiking, y.bs_mean_spiking1);
			data.mean_spiking_ds      = cat(1, data.mean_spiking_ds, y.bs_mean_spiking_ds1);
			data.flash_count          = cat(1, data.flash_count, y.bs_flash_count1);
			data.flash_count_ds       = cat(1, data.flash_count_ds, y.bs_flash_count_ds1);

			metadata.monkey           = cat(1, metadata.monkey, repmat(monkeys{m}, neurons, 1));
			metadata.session_date     = cat(1, metadata.session_date, repmat(y.session_date, neurons, 1));
			metadata.channel_num      = cat(1, metadata.channel_num, y.channel_num);
			metadata.iso_qual         = cat(1, metadata.iso_qual, y.iso_qual);
			metadata.depth            = cat(1, metadata.depth, y.depth);
			metadata.waveform_width   = cat(1, metadata.waveform_width, y.waveform_width);
			metadata.num_flashes      = cat(1, metadata.num_flashes, y.num_flashes);

			% Combine data and metadata
			results = [];
			results.data = data;
			results.metadata = metadata;

			% Save data after each session
			save([savedir 'result_data.mat'], 'results', '-v7.3');

		end
	end


end