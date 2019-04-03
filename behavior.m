datadir = './spatial_category_exp/';

colors = 3;
hlocs = 13;
vlocs = 26;
monkeys = {'quincy';'wahwah'};
 
all_lever_response = zeros(colors,hlocs,vlocs);
all_total_count = zeros(colors,hlocs,vlocs);
all_zone_color_response = zeros(colors,colors);
all_zone_color_count = zeros(colors,colors);

for m = 1:length(monkeys);
	
	if strcmp(monkeys{m}, 'wahwah')
		data_files = dir([datadir 'wahwah*Oct*.mat']);
	elseif strcmp(monkeys{m}, 'quincy')
		data_files = dir([datadir 'quincy*.mat']);
	end
	
	for i = 1:5; %length(data_files);
		load([datadir data_files(i).name])
		disp(['Session ' data_files(i).name ' loaded'])
		[lever_response, total_count, zone_color_response, zone_color_count] = calculate_response(session, monkeys{m});
		all_lever_response = all_lever_response + lever_response;
		all_total_count = all_total_count + total_count;
		all_zone_color_response = all_zone_color_response + zone_color_response;
		all_zone_color_count = all_zone_color_count + zone_color_count;

	end
end

lever_response = all_lever_response;			% 3 x 13 x 26
total_count = all_total_count;					% 3 x 13 x 26
zone_color_response = all_zone_color_response;	% 3 x 3
zone_color_count = all_zone_color_count;		% 3 x 3


%%%%% Lever Pulls at Each Location %%%%%

rp = squeeze(lever_response(1,:,:)./total_count(1,:,:));
gp = squeeze(lever_response(2,:,:)./total_count(2,:,:));
wp = squeeze(lever_response(3,:,:)./total_count(3,:,:));

data = {rp gp wp};
name = {'Red' 'Green' 'White'};

figure
for i=1:3
	curr_data = data{i};
	subplot(1,4,i),imagesc(0:12,-12.5:12.5,curr_data',[0,1]);
	axis equal;
	axis tight;
	colorbar;
	xlabel('X-location')
	ylabel('Y-location')
	title(strcat(name{i}, ' Flashes'))
end


%%%%% Color Responses Per Zone %%%%%

subplot(1,4,4),bar(zone_color_response./zone_color_count);
legend({'Red Flashes';'Green Flashes';'White Flashes'});
set(gca,'xticklabel',{'Red Zone';'Green Zone';'White Zone'},'ylim',[0,1]);
ylabel('Response Fraction');
title('Color Response');


% print('FillPageFigure','-dpdf','-fillpage')



function [lever_response, total_count, zone_color_response, zone_color_count] = calculate_response(sess, monkey)

	colors = 3;
	hlocs = 13;
	vlocs = 26;

	response_window = [150 400];
	num_trials = length(sess.x);

	lever_response = zeros(colors,hlocs,vlocs);
	total_count = zeros(colors,hlocs,vlocs);

	zone_color_response = zeros(3,3);
	zone_color_count = zeros(3,3);
	b = 1; % block number

	num_dead_flashes = 1;

	for i = 1:num_trials
		
		trial = sess.x(i);
		if strcmp(trial.trial_type,'ball_strikes') & trial.block_number == b & length(trial.flash_time) > 1

			trial.flash_pos(:,1) = trial.flash_pos(:,1) + 1;
			trial.flash_pos(:,2) = trial.flash_pos(:,2) + 13.5;

			% use flashes up to target flash, or second to last flash
			target_flash_ind = find((trial.flash_color==1 & trial.flash_zone==1) | (trial.flash_color==2) & (trial.flash_zone==2));
			if isempty(target_flash_ind)
				last_flash = length(trial.flash_time) - 1;
			else
				last_flash = target_flash_ind;
			end

			for j = num_dead_flashes+1:last_flash
				if isempty(trial.fixation_break) || trial.fixation_break > trial.flash_time(j) + response_window(2)
					c_ind = trial.flash_color(j);
					z_ind = trial.flash_zone(j);
					x_ind = trial.flash_pos(j,1);
					y_ind = trial.flash_pos(j,2);
					total_count(c_ind, x_ind, y_ind) = total_count(c_ind, x_ind, y_ind) + 1;
					zone_color_count(c_ind, z_ind) = zone_color_count(c_ind, z_ind) + 1;
					if ~isempty(trial.lever_up) ...
						&& trial.lever_up > trial.flash_time(j) + response_window(1) ...
						&& trial.lever_up <= trial.flash_time(j) + response_window(2)
						lever_response(c_ind, x_ind, y_ind) = lever_response(c_ind, x_ind, y_ind) + 1;
						zone_color_response(c_ind, z_ind) = zone_color_response(c_ind, z_ind) + 1;
					end
				end
			end
		end
	end
end

