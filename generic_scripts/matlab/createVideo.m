function createVideo(frame_array,frame_rate, x_axis, x_limit, y_axis, y_limit, dynamic_range, file_name, video_title)
%CREATE_VIDEO Summary of this function goes here
%   Detailed explanation goes here
% create the video writer
    
     fig = figure
     writerObj = VideoWriter(file_name);
     writerObj.FrameRate = frame_rate;      
     nFrames = size(frame_array,2);
     open(writerObj)   
        

     for i=1:size(frame_array,2)   
        i
        frame = 20*log10(abs(frame_array{i}./max(frame_array{i}(:))));

        fig = transpose(imagesc(x_axis,y_axis,frame, [-dynamic_range 0]));
        hold on
%         title(video_title + ' - CPI no.: ' + i)
        ylabel('Doppler Velocity [m/s]')
        xlabel('Range (m)')   
        c = colorbar;
        c.Label.String='Norm Power [dB]';

        ylim([-y_limit y_limit]); % Doppler Limit
        xlim([-inf x_limit]);   % Range Limit
        
        colormap jet;
        if ismac
            % Code to run on Mac platform
        elseif isunix
            saveas(fig,'/tmp/frame','png');
            im = imread('/tmp/frame.png');
        elseif ispc
            saveas(fig, tempdir + "frame",'png');
            im = imread(tempdir + "frame",'png');
        else
            disp('Platform not supported')
        end
        im = imresize(im,[875 1200]);
        writeVideo(writerObj, im);
     end
 
 close(writerObj);

end