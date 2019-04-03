function get_channel_nums()

	addpath('./Support Code/')

	datadir = './Processed Data/';
	savedir = './testing_result_data/';
	monkeys = {'quincy';'wahwah'};

	metadata = [];
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
			data_files = dir([datadir 'Wahwahball_strikes_*Oct*.mat']);
		elseif strcmp(monkeys{m}, 'quincy')
			data_files1 = dir([datadir 'Quincy_ball_strikes_*Feb*.mat']);
			data_files2 = dir([datadir 'Quincy_ball_strikes_*Mar*.mat']);
			data_files = [data_files1 ; data_files2];
		end
		
		% Iterate over the selected files
		for i = 1:length(data_files);
			
			% Load and process a session
			load([datadir data_files(i).name])
			disp(['Session ' data_files(i).name ' loaded'])

			y = [];
			y.session_date = session.ExperimentDate;
			y.channel_num = cat(1,session.NeuronInfo(:).Channel);
			y.iso_qual = cat(1,session.NeuronInfo(:).IsolationQuality);
			y.depth = cat(1,session.NeuronInfo(:).Depth);
			y.waveform_width = analyze_waveforms(session)';
			
			neurons = length(session.x(1).spike_times);
			for n=1:neurons
				y.num_flashes(n,:) = size(y.single_trial_combined(n).spike_resp, 1);
			end


			metadata.monkey           = cat(1, metadata.monkey, repmat(monkeys{m}, neurons, 1));
			metadata.session_date     = cat(1, metadata.session_date, repmat(y.session_date, neurons, 1));
			metadata.channel_num      = cat(1, metadata.channel_num, y.channel_num);
			metadata.iso_qual         = cat(1, metadata.iso_qual, y.iso_qual);
			metadata.depth            = cat(1, metadata.depth, y.depth);
			metadata.waveform_width   = cat(1, metadata.waveform_width, y.waveform_width);
			metadata.num_flashes      = cat(1, metadata.num_flashes, y.num_flashes);

			save([savedir 'metadata_old_code.mat'], 'metadata', '-v7.3');
		end
	end
end