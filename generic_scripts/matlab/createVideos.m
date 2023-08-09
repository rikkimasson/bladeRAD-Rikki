function createVideo(active_array,passive_array,frame_rate,active_range_axis, passive_range_axis, active_doppler_axis,passive_doppler_axis,range_limit,doppler_limit,dynamic_range,file_name,video_title)
%CREATE_VIDEO Summary of this function goes here
%   Detailed explanation goes here
% create the video writer
     writerObj = VideoWriter(file_name);
     writerObj.FrameRate = frame_rate;           
     open(writerObj)   
     fig = figure();
%     [A B C D] = fig.Position;
    fig.Position(3:4) = [1400 700];


%      fig.WindowState = 'maximized';

     for i=1:size(active_array,2)-1   

        active_frame =  10*log10(abs(active_array{i}./max(active_array{i}(:))));
        passive_frame = 10*log10(abs(passive_array{i}./max(passive_array{i}(:))));

        
        
        fig = subplot(1,2,1);
        imagesc(active_range_axis,active_doppler_axis, active_frame, [-dynamic_range, 0]);
        %title("Pre-CLEANed Passive Radar Capture")
        title("FMCW Active Radar Capture")
        ylabel('Doppler (Hz)')
        ylabel('Velocity (m/s)')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-doppler_limit doppler_limit]);
        xlim([-inf range_limit]);
%         colormap jet;

        
        subplot(1,2,2);
        imagesc(passive_range_axis, passive_doppler_axis, passive_frame, [-dynamic_range, 0]);
        %title('CLEANed Passive Radar Capture')
        title("Passive Radar Capture")
        ylabel('Doppler (Hz)')
        ylabel('Velocity (m/s)')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-doppler_limit doppler_limit]);
        xlim([-inf range_limit]);
%         colormap jet;
        
        
        saveas(fig,'/tmp/frame','png');
        im = imread('/tmp/frame.png');
        im = imresize(im,[1187 1094]);
        writeVideo(writerObj, im);
     end
 
 close(writerObj);

end