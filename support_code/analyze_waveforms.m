function [waveform_width, waveform_amp, waveform_half_width, waveform_half_width_neg] = analyze_waveforms(session)

waveform_width = [];
waveform_amp = [];
waveform_half_width = [];
waveform_half_width_neg = [];

try

snipet_size = size(session.NeuronInfo(1).Waveforms, 2);
min_cap = snipet_size*25-750; % only include spikes in which the min values occurs before this value


for j = 1:length(session.NeuronInfo)
    
    if size(session.NeuronInfo(j).Waveforms, 2) > 1
        [~, ix] = max(session.NeuronInfo(j).Waveforms(:,1:end)'); % discard the first eight time points, min or max should not occur within first ten points
        [~, iy] = min(session.NeuronInfo(j).Waveforms(:,1:end)');
        
        num_spikes = size(session.NeuronInfo(j).Waveforms, 1);
        aligned_wavforms = zeros(num_spikes, 876, 'single');
        count = 0;
        
        if num_spikes < 5000
            u = 1:num_spikes;
        else
            q = randperm(5000);
            u = q(1:5000);
        end
%         temp = [];
        for n = u
            s = interp(double(session.NeuronInfo(j).Waveforms(n,:)),25);
            [min_val, ind_min] = min(s(126:end));
            ind_min = ind_min + 125;
            [max_val, ind_max] = max(s(ind_min+1:end));
            ind_max = ind_max + ind_min;
            d1 = diff(s(1:end-1));
            d2 = diff(s(2:end));
%             d3 = diff(s(ind_min:ind_max),2);
            
%             num_inflection_pts = sum(d3(1:end-1).*d3(2:end)<0);
            
%             if (ind_min > 5*25 & ind_min < min_cap) && (all(d1(ind_min+1:ind_max-2).* d2(ind_min+1:ind_max-2)>0) & num_inflection_pts <= 1) 
            if (ind_min > 5*25 & ind_min < min_cap) && (all(d1(ind_min+1:ind_max-2).* d2(ind_min+1:ind_max-2)>0))
                count = count + 1;
                aligned_wavforms(count,:) = s(ind_min-125:ind_min+750);
%                 temp = cat(1,temp, ind_max-ind_min);
            end
        end

    else
        % for session in which only the template is saced, this will assume
        % two waveforms exists, both exact copies of the template. Bit of a
        % kludge
        [~, ix] = max(session.NeuronInfo(j).Waveforms(1:end)');
        [~, iy] = min(session.NeuronInfo(j).Waveforms(1:end)');
        session.NeuronInfo(j).Waveforms = session.NeuronInfo(j).Waveforms';
        aligned_wavforms = repmat(session.NeuronInfo(j).Waveforms(1:end),2,1);
        count = 1;
    end
    
    
    aligned_wavforms = aligned_wavforms(1:count, :);
    y = mean(aligned_wavforms);
  
    [max_val, ind_max] = max(y);
    [min_val, ind_min] = min(y);
    m = max_val/2;
    
    [~,ix] = min(abs(y(ind_min+1:ind_max)-m));
    [~,iy] = min(abs(y(1+ind_max:end)-m));
    ix = ix + ind_min;
    iy = iy + ind_max;
    
    waveform_width = [waveform_width; ind_max-ind_min];
    waveform_amp = [waveform_amp; max_val];
    if ~isempty(iy) & ~isempty(ix)
        waveform_half_width = [waveform_half_width; iy-ix];
    else
        waveform_half_width = [waveform_half_width; -1];
    end
    
    
    m = min_val/2;
    [~,ix] = min(abs(y(1:ind_min)-m));
    [~,iy] = min(abs(y(1+ind_min:end)-m));
    iy = iy + ind_min;
    
    if ~isempty(iy) & ~isempty(ix)
        waveform_half_width_neg = [waveform_half_width_neg; iy-ix];
    else
        waveform_half_width_neg = [waveform_half_width_neg; -1];
    end   
end

catch
end