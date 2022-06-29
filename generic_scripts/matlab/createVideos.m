function createVideo(active_array,passive_array,frame_rate,active_range_axis, passive_range_axis, active_doppler_axis,passive_doppler_axis,dynamic_range,file_name,video_title)
%CREATE_VIDEO Summary of this function goes here
%   Detailed explanation goes here
% create the video writer
     writerObj = VideoWriter(file_name);
     writerObj.FrameRate = frame_rate;           
     open(writerObj)   
     fig = figure();
     fig.WindowState = 'maximized';

     for i=1:size(active_array,2)-1   

        active_frame =  10*log(abs(active_array{i}./max(active_array{i}(:))));
        passive_frame = 10*log(abs(passive_array{i}./max(passive_array{i}(:))));

        
        
        fig = subplot(1,2,1);
        imagesc(active_range_axis,active_doppler_axis, active_frame, [-dynamic_range, 0]);
        %title("Pre-CLEANed Passive Radar Capture")
        title("FMCW Active Radar Capture")
        ylabel('Doppler (Hz)')
        ylabel('Velocity (m/s)')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-100 100]);
        ylim([-10 10]);
        xlim([-inf 50]);
        colormap jet;

        
        subplot(1,2,2);
        imagesc(passive_range_axis, passive_doppler_axis, passive_frame, [-dynamic_range, 0]);
        %title('CLEANed Passive Radar Capture')
        title("Wi-Fi Passive Radar Capture")
        ylabel('Doppler (Hz)')
        ylabel('Velocity (m/s)')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power (dB)';
        ylim([-100 100]);
        ylim([-10 10]);
        xlim([-inf 50]);
        colormap jet;
        
        
        saveas(fig,'frame','png');
        im = imread('frame.png');
        writeVideo(writerObj, im);
     end
 
 close(writerObj);

end