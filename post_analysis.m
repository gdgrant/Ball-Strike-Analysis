function post_analysis()

	addpath('./support_code/')
	load('./result_data/result_data.mat')

	neuron_pairs = pair_nearby_neurons(results.metadata);


end