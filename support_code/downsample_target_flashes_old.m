function [spike_resp_ds, flash_col_ds, flash_pos_ds, flash_zone_ds, target_flash_ds] = downsample_target_flashes_old(spike_resp, flash_pos, flash_col, flash_zone, green_or_cyan)

spike_resp_ds = [];
flash_pos_ds = [];
flash_col_ds = [];
target_flash_ds = [];
flash_zone_ds = [];

% for this calculation, we just want to ensure that we're not
% mixing target and non-target flahses
% flash_zone(flash_col == 1 & flash_zone~=1) = 0;
% flash_zone(flash_col == 2 & flash_zone~=2) = 0;

% for white flahses, we don't care about flash zone
% flash_zone(flash_col == 3) = 0;

flash_zone(flash_zone>=3) = 3;

% use all flashes not in the target zone and calculate the flash count
% histogram
for col = 1:3
    falsh_count_hist = [];
    for x = 1:13
        for y = 1:26
            current_zone = get_zone(x-1, y-13.5, green_or_cyan);
            % check whether it's a target location
            if current_zone ~= col | col > 2
                ind = find(flash_col == col & flash_pos(:,1)==x & flash_pos(:,2)==y);
                falsh_count_hist = cat(1, falsh_count_hist, length(ind));
            end
        end
    end
    
    
%     target_flashes_to_use_index = randperm(length(falsh_count_hist));
%     target_flash_num = 0;
    % select only a subset of  flahses based on the flash count histogram
    % calculated for non-target flahses
    
    used_ind = []; % will keep track of which flashes are used so we won't have any duplicates
    
    % select grid points at random
    [x_grid,y_grid]= meshgrid(1:13,1:26);
    num_grid_points = numel(x_grid);
    x_grid = reshape(x_grid,num_grid_points,1);
    y_grid = reshape(y_grid,num_grid_points,1);
    q = randperm(num_grid_points);
    x_grid = x_grid(q);
    y_grid = y_grid(q);
    
    for i = 1:num_grid_points
        x = x_grid(i);
        y = y_grid(i);
        
        current_zone = get_zone(x-1, y-13.5, green_or_cyan);
%         current_zone = min(3, current_zone);
        ind = find(flash_col == col & flash_pos(:,1)==x & flash_pos(:,2)==y);
        ind = setdiff(ind, used_ind);
        %             if ~isempty(ind)
        
        % by using this if statement, we're only downsampling target flahses;
        % downsampling non-target flahses leads to a negative bias
        % inside the target zones
        %                 if current_zone==col & col < 3
        % select number of target flahses with replacement
        q = randperm(length(falsh_count_hist));
        flahses_to_use = falsh_count_hist(q(1));
        % select number of target flahses without replacement
        %                 target_flash_num = target_flash_num + 1;
        %                 flahses_to_use = falsh_count_hist(target_flashes_to_use_index(target_flash_num));
        %                 else
        %                     flahses_to_use = length(ind);
        %                 end
        
        if length(ind) >= flahses_to_use
            % random without replacement
            q = randperm(length(ind));
            ind = ind(q(1:flahses_to_use));
            % random with replacement
            %                         q = randi(length(ind),1,flahses_to_use);
            %                         ind = ind(q);
        else
            
            % if there are not enough flashes, use flashes from
            % adjacent points
            num_additional_flashes = flahses_to_use - length(ind);
            %                     ind2 = find(flash_col == col & flash_pos(:,1)>=x-1 & flash_pos(:,1)<=x+1 & flash_pos(:,2)>=y-1 ...
            %                         & flash_pos(:,2)<=y+1 & flash_zone==current_zone);
            ind2 = find(flash_col == col & flash_pos(:,1)>=x-1 & flash_pos(:,1)<=x+1 & flash_pos(:,2)>=y-1 ...
                & flash_pos(:,2)<=y+1);
            ind2 = setdiff(ind2, used_ind);
            ind2 = setdiff(ind2, ind);
            
            
            
            if length(ind2)>=num_additional_flashes
                %                             % random without replacement
                q = randperm(length(ind2));
                ind2 = ind2(q(1:num_additional_flashes));
                %                             % random with replacement
                % %                             q = randi(length(ind2),1,num_additional_flashes);
                % %                             ind2 = ind2(q);
            else
                % if still not enough flashes, extend search radius
                %                         ind2 = find(flash_col == col & flash_pos(:,1)>=x-2 & flash_pos(:,1)<=x+2 & flash_pos(:,2)>=y-2 ...
                %                             & flash_pos(:,2)<=y+2 & flash_zone==current_zone);
                ind2 = find(flash_col == col & flash_pos(:,1)>=x-2 & flash_pos(:,1)<=x+2 & flash_pos(:,2)>=y-2 ...
                    & flash_pos(:,2)<=y+2);
                
                ind2 = setdiff(ind2, used_ind);
                ind2 = setdiff(ind2, ind);
                q = randperm(length(ind2));
                num_additional_flashes = min(num_additional_flashes, length(ind2));
                ind2 = ind2(q(1:num_additional_flashes));
                
            end
            ind = union(ind, ind2);
        end
        if current_zone==col & col < 3
            is_target_flash = 1;
        else
            is_target_flash = 0;
        end
        
        
        % add flahses to list of flashes we've added
        if ~isempty(ind)
            if size(ind,2)>1
                ind = ind';
            end
            used_ind = cat(1,used_ind, ind);
        end
        
        flash_pos_ds = cat(1, flash_pos_ds, repmat([x y],length(ind),1));
        flash_col_ds = cat(1, flash_col_ds, repmat(col,length(ind),1));
        target_flash_ds = cat(1, target_flash_ds, is_target_flash*ones(length(ind),1));
        flash_zone_ds = cat(1, flash_zone_ds, repmat(current_zone,length(ind),1));
        spike_resp_ds = cat(1, spike_resp_ds, spike_resp(ind,:));
        
        %             end
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
