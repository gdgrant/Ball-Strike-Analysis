function waveform_width = cull_waveforms(session)

waveform_width = [];

for j = 1:length(session.NeuronInfo)
    if size(session.NeuronInfo(j).Waveforms, 2) > 1
        [~, ix] = max(session.NeuronInfo(j).Waveforms(:,11:end)'); % discard the first ten time points, min or max should not occur within first ten points
        [~, iy] = min(session.NeuronInfo(j).Waveforms(:,11:end)');
    else
        [~, ix] = max(session.NeuronInfo(j).Waveforms(11:end)');
        [~, iy] = min(session.NeuronInfo(j).Waveforms(11:end)');
    end
    d = ix - iy;
    md = median(d);
    good_spike_ind = find(abs(d - md) <= 4); % only count those spikes for which the width is within 4 steps of the median
    mean_width = mean(d(good_spike_ind));
    waveform_width = [waveform_width mean_width];
end