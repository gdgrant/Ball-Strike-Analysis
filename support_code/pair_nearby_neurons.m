function [neuron_pairs, channel_nums] = pair_nearby_neurons(metadata)

	monkey   = metadata.monkey;
	channel  = metadata.channel_num;
	iso_qual = metadata.iso_qual;
	depth    = metadata.depth;

	monkeys = {"quincy", "wahwah"};
	iso_qual_mask = iso_qual >= 3;

	neuron_pairs = [];	% [2 (inds) x pairs]
	channel_nums = [];	% [2 (monkey) x pairs]


	for m=1:2

		monkey_mask = strcmp(monkey, monkeys{m});

		for c=1:max(channel)

			channel_mask = channel==c;

			neuron_inds = find(iso_qual_mask & monkey_mask & channel_mask);
			num_inds = length(neuron_inds);

			% If there are not enough neurons in the channel to make a pair,
			% continue to the next channel
			if num_inds < 2
				continue
			end

			% Iterate over each possible neuron pair in the collected indices
			disp('-------')
			for i=1:num_inds
				for j=i:num_inds
					ind_i = neuron_inds(i);
					ind_j = neuron_inds(j);


					disp([i, j, depth(ind_i), depth(ind_j)])
				end
			end





		end
	end
end