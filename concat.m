function [data, metadata] = concat()

	files = [];
	files(1).name = 'result_data_m1_s23.mat';
	files(2).name = 'result_data_m2_s15.mat';

	% files = dir('./result_data/result_data_m*_s*.mat');

	data = [];
	metadata = [];

	data.perc_mask = [];
	data.time_windows = [];
	data.percentiles = [];

	metadata.monkey = [];
	metadata.session_date = [];
	metadata.channel_num = [];
	metadata.iso_qual = [];
	metadata.depth = [];
	metadata.waveform_width = [];

	for i = 1:length(files)
		load(['./result_data/' files(i).name])
		disp(['./result_data/' files(i).name])
		neurons = size(results.data.perc_mask, 1)

		data.perc_mask    = cat(1, data.perc_mask, single(results.data.perc_mask));
		data.time_windows = cat(1, data.time_windows, results.data.time_windows);
		data.percentiles  = cat(1, data.percentiles, results.data.percentiles);

		metadata.monkey         = cat(1, metadata.monkey, results.metadata.monkey);
		metadata.session_date   = cat(1, metadata.session_date, results.metadata.session_date);
		metadata.channel_num    = cat(1, metadata.channel_num, results.metadata.channel_num);
		metadata.iso_qual       = cat(1, metadata.iso_qual, results.metadata.iso_qual);
		metadata.depth          = cat(1, metadata.depth, results.metadata.depth);
		metadata.waveform_width = cat(1, metadata.waveform_width, results.metadata.waveform_width);

		whos
	end

end