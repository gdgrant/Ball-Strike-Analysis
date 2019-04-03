function [flash_lever_response, flash_count, flash_zone_color_response, flash_zone_color_count, z] = ball_strikes_behavior


addpath('C:\Users\nicol_000\Projects\MatlabTools\Error Bar Plots')
addpath('C:\Users\nicol_000\Projects\Ball Strikes Experiment\Code\Support Code')
data_dir = 'C:\Users\nicol_000\Projects\Ball Strikes Experiment\Processed Data\';


flash_lever_response = zeros(3, 13, 26);
flash_count = zeros(3,13,26);
flash_zone_color_response = [];
flash_zone_color_count = [];

monkey = {'quincy';'wahwah'};

for m = 1:2
    if strcmp(monkey{m}, 'wahwah')
        data_files = dir([data_dir 'Wahwahball_strikes_*Oct*reduced.mat']);
    elseif strcmp(monkey{m}, 'quincy')
        data_files = dir([data_dir 'Quincy_ball_strikes_*reduced.mat']);
    end
    
    z = [];
    
    for i = 1:length(data_files)
        load([data_dir data_files(i).name])
        disp(['Session ' data_files(i).name ' loaded'])
        [flr, fc, fzcr, fzcc] = calculate_distractor_response(session, monkey{m});
        flash_lever_response = flash_lever_response + flr;
        flash_count = flash_count + fc;
        
        flash_zone_color_response = cat(3,flash_zone_color_response,fzcr);
        flash_zone_color_count = cat(3,flash_zone_color_count,fzcc);
        
        x1=(fzcr(1,1)+fzcr(2,2))/(fzcc(1,1)+fzcc(2,2));
        x2=(fzcr(2,1)+fzcr(1,2))/(fzcc(2,1)+fzcc(1,2));
        z = [z; x1 x2];
    end 
end

% plot response rate bar plot
m = squeeze(mean(flash_zone_color_response./flash_zone_color_count,3));
se = squeeze(std(flash_zone_color_response./flash_zone_color_count,[],3))/sqrt(size(flash_zone_color_count,3));
errorbar_groups(100*m,100*se)
set(gca,'xticklabel',{'Red Flash';'Green Flash';'White Flash'})
ylabel('Response rate (%)')
legend({'Flash inside RED zone';'Flash inside GREEN zone';'Flash outside both target zones'})


% create response heat maps
flash_titles = {'Red Flashes';'Green Flashes';'White Flashes'};
figure
for j = 1:3
    subplot(1,3,j)
    imagesc(0:12,-12.5:12.5,squeeze(100*flash_lever_response(j,:,:)./flash_count(j,:,:))',[0 100])
    axis xy
    colorbar
    title(flash_titles{j})
    xlabel('X-axis (deg)')
    ylabel('Y-axis (deg)')
    colormap(cool)
    hold on
    y1 = -0.1;y2 = 0.1;
    plot([4 8],[y1 y1],'g','LineWidth',2);plot([4 8],[-5 -5],'g','LineWidth',2);plot([4 4],[-5 y1],'g','LineWidth',2);plot([8 8],[-5 y1],'g','LineWidth',2)
    plot([4 8],[y2 y2],'r','LineWidth',2);plot([4 8],[5 5],'r','LineWidth',2);plot([4 4],[5 y2],'r','LineWidth',2);plot([8 8],[5 y2],'r','LineWidth',2)
end


% create projection plots
figure
subplot(1,2,1)
plot(-12.5:12.5,100*squeeze(mean(flash_lever_response(1,4:8,:)./flash_count(1,4:8,:))),'r','LineWidth',2)
hold on
plot(-12.5:12.5,100*squeeze(mean(flash_lever_response(2,4:8,:)./flash_count(2,4:8,:))),'g','LineWidth',2)
plot([0.1 0.1],[0 100],'r--','LineWidth',1)
plot([5 5],[0 100],'r--','LineWidth',1)
plot([-0.1 -0.1],[0 100],'g--','LineWidth',1)
plot([-5 -5],[0 100],'g--','LineWidth',1)
axis([-12.5 12.5 0 100])
xlabel('Y-axis (deg)')
ylabel('Response rate (%)')

subplot(1,2,2)
plot(0:12,100*squeeze(mean(flash_lever_response(1,:,9:16)./flash_count(1,:,9:16),3)),'r','LineWidth',2)
hold on
plot(0:12,100*squeeze(mean(flash_lever_response(2,:,9:16)./flash_count(2,:,9:16),3)),'g','LineWidth',2)
plot([3.45 3.45],[0 100],'r--','LineWidth',1)
plot([8.45 8.45],[0 100],'r--','LineWidth',1)
plot([3.55 3.55],[0 100],'g--','LineWidth',1)
plot([8.55 8.55],[0 100],'g--','LineWidth',1)
axis([0 12 0 100])
xlabel('X-axis (deg)')
ylabel('Response rate (%)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_behavior(flash_lever_response, flash_count, flash_zone_color_response, flash_zone_color_count)

r1=sum(flash_lever_response(1,:,14:18),3)./sum(flash_count(1,:,14:18),3);
g1=sum(flash_lever_response(2,:,9:13),3)./sum(flash_count(2,:,9:13),3);
r2=squeeze(sum(flash_lever_response(1,4:8,:),2)./sum(flash_count(1,4:8,:),2));
g2=squeeze(sum(flash_lever_response(2,4:8,:),2)./sum(flash_count(2,4:8,:),2));

figure;subplot(1,2,1),
plot(r1,'-ro')
hold on
plot(g1,'-go')
xlabel('X-axis location')
ylabel('Response fraction')
subplot(1,2,2),
plot(-12.5:12.5,r2,'-ro')
hold on
plot(-12.5:12.5,g2,'-go')
xlabel('Y-axis location')
ylabel('Response fraction')
subplot(1,2,1),
 plot([3.5 3.5],[0 1],'k--')
plot([8.5 8.5],[0 1],'k--')
subplot(1,2,2),
plot([0 0],[0 1],'k--')
plot([-5 -5],[0 1],'k--')
plot([5 5],[0 1],'k--')
% suptitle('Quincy Response Rate to Red and Green Flashes')


figure
bar(flash_zone_color_response./flash_zone_color_count)
legend({'Red Flashes';'Green Flashes';'White Flashes'})
set(gca,'xticklabel',{'Red Zone';'Green Zone';'White Zone'})
ylabel('Response fraction')
% title('Quincy Response Rate')

r1=squeeze(flash_lever_response(1,:,:)./flash_count(1,:,:));
g1=squeeze(flash_lever_response(2,:,:)./flash_count(2,:,:));
w1=squeeze(flash_lever_response(3,:,:)./flash_count(3,:,:));

figure
subplot(1,3,1),imagesc(0:12,-12.5:12.5,r1',[0 1]);colorbar;xlabel('X-location');ylabel('Y-location');title('Red Flashes')
subplot(1,3,2),imagesc(0:12,-12.5:12.5,g1',[0 1]);colorbar;xlabel('X-location');ylabel('Y-location');title('Green Flashes')
subplot(1,3,3),imagesc(0:12,-12.5:12.5,w1',[0 1]);colorbar;xlabel('X-location');ylabel('Y-location');title('White Flashes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flash_lever_response, flash_count, flash_zone_color_response, flash_zone_color_count] = calculate_distractor_response(session, monkey)

response_window = [150 400];
num_trials = length(session.x);

flash_lever_response = zeros(3, 13, 26);
flash_count = zeros(3,13,26);

flash_zone_color_response = zeros(3,3);
flash_zone_color_count = zeros(3,3);
b = 1; % desired block number

if strcmp(monkey, 'wahwah')
    num_dead_flashes = 1;
elseif strcmp(monkey, 'quincy')
    num_dead_flashes = 1;
end



for i = 1:num_trials
    if strcmp(session.x(i).trial_type,'ball_strikes') & session.x(i).block_number == b
        
        [flash_time, flash_color, flash_zone, flash_pos] = convert_flash_sequence(session.x(i), 'green');
        if length(flash_time) > 1
            flash_pos(:,1) = flash_pos(:,1) + 1;
            flash_pos(:,2) = flash_pos(:,2) + 13.5;
            
            % only use flashes up until the target flash, or the second to
            % last flash
            target_flash_ind = find((flash_color==1 & flash_zone==1) | (flash_color==2 & flash_zone==2));
            if ~isempty(target_flash_ind)
                last_flash = target_flash_ind;
            else
                last_flash = length(flash_time) - 1;
            end
            
            for j = num_dead_flashes+1:last_flash
                if isempty(session.x(i).fixation_break) || session.x(i).fixation_break > flash_time(j) + response_window(2)
                    flash_count(flash_color(j), flash_pos(j,1), flash_pos(j,2)) = flash_count(flash_color(j), flash_pos(j,1), flash_pos(j,2)) + 1;
                    flash_zone_color_count(flash_color(j),flash_zone(j)) = flash_zone_color_count(flash_color(j),flash_zone(j)) + 1;
                    if ~isempty(session.x(i).lever_up) && session.x(i).lever_up > flash_time(j) + response_window(1) & session.x(i).lever_up <= flash_time(j) + response_window(2)
                        flash_lever_response(flash_color(j), flash_pos(j,1), flash_pos(j,2)) = flash_lever_response(flash_color(j), flash_pos(j,1), flash_pos(j,2)) + 1;
                        flash_zone_color_response(flash_color(j),flash_zone(j)) = flash_zone_color_response(flash_color(j),flash_zone(j)) + 1;
                    end
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [flash_time, flash_color, flash_zone, flash_pos] = convert_flash_sequence(x)


% flash_time = [x.red_distractor_flash; x.green_distractor_flash; x.white_distractor_flash];
% flash_pos = [x.red_flash_pos; x.green_flash_pos; x.white_flash_pos];
% flash_color = [ones(length(x.red_distractor_flash), 1); 2*ones(length(x.green_distractor_flash), 1); 3*ones(length(x.white_distractor_flash), 1)];
% if ~isempty(x.target_flash)
%     flash_time = [flash_time; x.target_flash];
%     if ~isempty(x.green_target_pos)
%         flash_pos = [flash_pos; x.green_target_pos];
%     elseif ~isempty(x.red_target_pos)
%         flash_pos = [flash_pos; x.red_target_pos];
%     else
%         error('Something went wrong.')
%     end
%     target_zone = get_zone(flash_pos(end,:));
%     if target_zone == 1
%         flash_color = [flash_color; 1];
%     elseif target_zone == 2
%         flash_color = [flash_color; 2];
%     else
% %         keyboard
% %         error('Something went wrong.')
%     end
% end
% flash_zone = [];
% if size(flash_pos, 1) > 0
%     for i = 1:length(flash_time)
%         zone = get_zone(flash_pos(i,:));
%         flash_zone = [flash_zone; zone];
%     end
%     [~,ind] = sort(flash_time, 'ascend');
%     flash_time = flash_time(ind);
%     flash_zone = flash_zone(ind);
%     flash_color = flash_color(ind);
%     flash_pos = flash_pos(ind, :);
% else
%     flash_time = [];
%     flash_zone = [];
%     flash_color = [];
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function zone = get_zone(flash_location)
% 
% green_targetX = [4 8];
% green_targetY = [-4.5 -0.5];
% red_targetX = [4 8];
% red_targetY = [0.5 4.5];
% 
% if flash_location(1) >= green_targetX(1) & flash_location(1) <= green_targetX(2) ...
%         & flash_location(2) >= green_targetY(1) & flash_location(2) <= green_targetY(2)
%     zone = 2;
% elseif flash_location(1) >= red_targetX(1) & flash_location(1) <= red_targetX(2) ...
%         & flash_location(2) >= red_targetY(1) & flash_location(2) <= red_targetY(2)
%     zone = 1;
% else
%     zone = 3;
% end
