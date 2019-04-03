function session = demo_response(session)


	trials = length(session.x);
	neurons = length(session.x(1).spike_times);

	for t=1:trials
		flashes = length(session.x(t).flash_time);

		if strcmp(session.x(t).trial_type, 'ball_strikes') & (flashes > 0)
			for f=1:flashes

				time = session.x(t).flash_time(f);
				xpos = session.x(t).flash_pos(f,1);
				ypos = session.x(t).flash_pos(f,2);
				col = session.x(t).flash_color(f);

				if xpos == 1.0 & ypos == 2.5 & col == 1
					for n=1:neurons
						session.x(t).spike_times{n} = cat(1,session.x(t).spike_times{n},time-100,time-50,time-25,time,time+25,time+50,time+100);
					end
				elseif xpos == 4.0 & ypos == 6.5 & col == 2
					for n=1:neurons
						session.x(t).spike_times{n} = cat(1,session.x(t).spike_times{n},time-100,time-50,time-25,time,time+25,time+50,time+100);
					end
				end

			end
		end
	end
end