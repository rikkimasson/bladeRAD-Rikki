function createCFARVideo(frame_array,frame_rate,x_axis,x_limit,y_axis,y_limit,file_name,video_title)
%CREATE_VIDEO Summary of this function goes here
%   Detailed explanation goes here
% create the video writer
     writerObj = VideoWriter(file_name);
     writerObj.FrameRate = frame_rate;           
     open(writerObj)   
     
     for i=1:size(frame_array,2)-1   

        fig =  transpose(imagesc(x_axis,y_axis, frame_array{i}));
        title(video_title +  " CPI Number: " + i)
        ylabel('Doppler (Hz)')
        xlabel('Range (m)')   
        ylim([-y_limit y_limit]);
        xlim([-inf x_limit]);
        saveas(fig,'frame','png');
        im = imread('frame.png');
        writeVideo(writerObj, im);
     end
 
 close(writerObj);

end