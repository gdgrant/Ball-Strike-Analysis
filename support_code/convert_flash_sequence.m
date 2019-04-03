function [flash_time, flash_color, flash_zone, flash_pos, flash_color_name, target_flash_time] = convert_flash_sequence(x, green_or_cyan)

flash_time = [x.red_distractor_flash; x.green_distractor_flash; x.white_distractor_flash];
flash_pos = [x.red_flash_pos; x.green_flash_pos; x.white_flash_pos];
flash_color = [ones(length(x.red_distractor_flash), 1); 2*ones(length(x.green_distractor_flash), 1); 3*ones(length(x.white_distractor_flash), 1)];
if ~isempty(x.target_flash)
    flash_time = [flash_time; x.target_flash];
    target_flash_time = x.target_flash;
    if ~isempty(x.green_target_pos)
        flash_pos = [flash_pos; x.green_target_pos];
    elseif ~isempty(x.red_target_pos)
        flash_pos = [flash_pos; x.red_target_pos];
    else
        error('Something went wrong.')
    end
    target_zone = get_zone(flash_pos(end,1), flash_pos(end,2), green_or_cyan);
    if target_zone == 1
        flash_color = [flash_color; 1];
    elseif target_zone == 2
        flash_color = [flash_color; 2];
    else
        keyboard
        error('Something went wrong.')
    end
else
    target_flash_time = [];
end
flash_zone = [];
if size(flash_pos, 1) > 0 && size(flash_pos, 1) == size(flash_time, 1) 
    for i = 1:length(flash_time)
%         zone = get_zone(flash_pos(i,:),green_or_cyan);
%         zone = get_zone_precise(flash_pos(i,:),green_or_cyan);
        zone = get_zone(flash_pos(i,1), flash_pos(i,2), green_or_cyan);
        flash_zone = [flash_zone; zone];
    end
    [~,ind] = sort(flash_time, 'ascend');
    flash_time = flash_time(ind);
    flash_zone = flash_zone(ind);
    flash_color = flash_color(ind);
    flash_pos = flash_pos(ind, :);
    flash_color_name = cell(1, length(flash_color));
    for i = 1:length(flash_color)
        if flash_color(i) == 1
            flash_color_name{i} = 'red';
        elseif flash_color(i) == 2
            flash_color_name{i} = 'green';
        elseif flash_color(i) == 3
            flash_color_name{i} = 'white';
        end
    end
else
    flash_time = [];
    flash_zone = [];
    flash_color = [];
    flash_color_name = cell(1, length(flash_color));
end

end

function zone = get_zone(x, y, green_or_cyan)

	% RED Target   - Zone 1
	% GREEN Target - Zone 2
	% OUTSIDE      - Zone 3

	if nargin < 3
		green_or_cyan = 'green';
	end

	red_zone_x = [4 8];
	red_zone_y = [0.5 4.5];

	if strcmp(green_or_cyan,'green')
		green_zone_x = [4 8];
		green_zone_y = [-4.5 -0.5];
	elseif strcmp(green_or_cyan,'cyan')
		green_zone_x = [6 10];
		green_zone_y = [5.5 9.5];
	end

	red_x = (red_zone_x(1) <= x) & (x <= red_zone_x(2));
	red_y = (red_zone_y(1) <= y) & (y <= red_zone_y(2));
	green_x = (green_zone_x(1) <= x) & (x <= green_zone_x(2));
	green_y = (green_zone_y(1) <= y) & (y <= green_zone_y(2));

	if red_x && red_y
		zone = 1;
	elseif green_x && green_y
		zone = 2;
	else
		zone = 3;
	end
end
